%% Spectroradiometer (PR670) Automatic measure code
% This code aims to test stability and additivity of the projector
% It basically measures lights automatically and continuously with a
% specific length of time interval

%% Initialize and Set paths
clear all; close all; clc;
tbUse('BrainardLabBase')
addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx')); % add VPixx toolbox ('Datapixx') to path (incl. all subfolders)

%% Start 1 (PR-670)
% In terminal, command to authorize to read and write 'sudo chmod a+rw /dev/ttyACM0 (the last character is number zero)'
% Mostly PR670 connects to 'ttyACM0', but sometimes to 'ttyACM1', so just
% do both to make it sure to be connected
% Also command to check the status 'ls -l /dev/ttyACM0'

command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
command_check = 'ls -l /dev/ttyACM0';
unix(command_check)

%% Start 2 (Vpixx projector)
% Run terminal command 'vputil' on here using unix function
% rw 0x1c8 0x1 % (Primary 0 always ON) > Screen goes red 
% rw 0x1c8 0x2 % (Primary 1 always ON) > Screen goes green
% rw 0x1c8 0x4 % (Primary 2 always ON) > Screen goes blue
% rw 0x1c8 0x7 % (All primaries on) > Screen goes brighter with no chromaticity changes (at least visually)
% rw 0x1c8 0x0 % (Set to default) > Back to the first screen

command_primaries = 'vputil rw 0x1c8 0x7 -q quit'; % Set to 'All Primaries on'
unix(command_primaries)

% Connect the Vpixx device
isReady = Datapixx('open');
isReady = Datapixx('IsReady'); % Returns non-0 if a Datapixx has been successfully opened for use.

% 'primaryColor' should be a value of either 0, 1, or 2
% 'subColor' is value between 0 and 15. / *Note: selecting subColor 8 does
% not have any effect (Unavailable)

%% Vpixx projector color control

% Set PrimaryColor(0,1,2) and SubColor(0-15 except 8)
current = 252; % 0-252

for j=1:3 % PrimaryColor
   for i=1:16 % SubColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current); 
   end
end

% Check the current settings
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

currents

%% *********ROUTINE STARTS FROM HERE********* (Stability)
% Measure the spectrum automatically in a specific time interval

% Set the measurement date (it reflects the saved file name)
Date = '0723'; 

% Time delay before start (go out the room before it's done)
timesec = 15;
t = timer('TimerFcn','stat=false; disp(''Time OFF to leave!'')','StartDelay',timesec);
start(t)
stat=true;
while(stat==true)
     disp('.')
     pause(1) % Pause for 1 seconds (unit)
end
delete(t)

% Connect PR670
CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201]; % Wavelength range (380-780 nm) with 2 nm interval

% Set the measurement time length in min
time_min = 60; % time in minutes
timedelay_min = 0.5;

time_sec = time_min*60;
timedelay_sec = timedelay_min*60;

NumMeasurements = time_min/timedelay_min + 1; % The number of measurement repetitions / +1 for the measurement at cold (right after turned on)

for i=1:NumMeasurements
    
    % PR670 measure
    fw(:,i)=MeasSpd(S,5,'all'); % Measurement (it takes a while)

    % Time delay
    t = timer('TimerFcn','stat=false; disp(''Measurement starts!'')','StartDelay',timedelay_sec);
    start(t)
    stat=true;
    while(stat==true)
    disp('.')
    pause(1) % Pause for 1 seconds (unit)
    end
    delete(t)
    
end

% Save the spectrum and XYZ values
filename = append(Date,'Spectrum','.mat'); 
save(filename,'fw') % Save the variable as a '.mat' file

% Save the XYZ values
% Wavelength linear interpolation (1/2/5 nm interval)
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range
fw_1nm = interp1(w_2nm,fw,w_1nm);
fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation
% xyY = XYZToxyY(XYZ'); % xyY
% colorgamut=XYZToxyY(T_XYZ');
% colorgamut(:,82)=colorgamut(:,1);

filename2 = append(Date,'XYZ','.mat'); 
save(filename2,'XYZ') 

%% *********ROUTINE STARTS FROM HERE********* (Additivity)
%% 1) Single Channel measurement (A total of 15 LED channels)
CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201]; % Wavelength range (380-780 nm) with 2 nm interval
Date = '0723';

for i=1:5 % Num of 15 LED channels

    % Set all channel currents to 0 as baseline
    isReady = Datapixx('open');
    isReady = Datapixx('IsReady');
    current_base = 0;
    for j=1:3 % PrimaryColor
        for k=1:16 % subColor
         Datapixx('SetPropixxHSLedCurrent', j-1, k-1, current_base);
        end
    end

    % Measurement with different currents
    % current_set = [10, 126, 252]; % current=10 will be used as 'black'
      current_set = [10:10:240,252]; % p.s. These all sets took 1hr 45min (a total of 368 measurements)

      for f = 1:length(current_set);
        current = current_set(f);

        % Measure each LED channel
        isReady = Datapixx('open');
        isReady = Datapixx('IsReady');
        % current = 252;
        Datapixx('SetPropixxHSLedCurrent', 0, i-1, current); % Primary 1
        Datapixx('SetPropixxHSLedCurrent', 1, i-1, current); % Primary 2
        Datapixx('SetPropixxHSLedCurrent', 2, i-1, current); % Primary 3

        % PR670 measurement
        fw_single(:,f+length(current_set)*(i-1))=MeasSpd(S,5,'all');
        disp(append('LED Channel ',num2str(i),' current ',num2str(f),' complete!'))
      end

end

% Save the spectrums (.mat file)
filename = append(Date,'Spectrum_single','.mat');
save(filename,'fw_single')% Wavelength linear interpolation (1/2/5 nm interval)

%% Stability (XYZ analysis)
% Save the XYZ values
% Wavelength linear interpolation (1/2/5 nm interval)
load('0721Spectrum(Stability).mat') % Load the data

t = [0:0.5:110.5]; % Time range in min
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range
fw_1nm = interp1(w_2nm,fw,w_1nm);
fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation

time_min = 20; % Set the time range for comparison
time_index = time_min*2;

figure(1); subplot(2,1,1); hold on;
plot(w_2nm,fw(:,1:time_index),'k');
plot(w_2nm,fw(:,time_index+1:222),'r');
xlabel('Wavelength (nm)');
ylabel('Spectral irradiance');
% legend('');

figure(1); subplot(2,1,2); hold on;
plot(t(1:time_index),XYZ(1:time_index,2),'k.-')
plot(t(time_index+1:222),XYZ(time_index+1:222,2),'r.-')
xlabel('Time (min)');
ylabel('Luminance (cd/m2)');


%% 2) Random combinations of the LED channels

Date = '0723';

% Time delay before start (go out the room before it's done)
timesec = 15;
t = timer('TimerFcn','stat=false; disp(''Time OFF to leave!'')','StartDelay',timesec);
start(t)
stat=true;
while(stat==true)
     disp('.')
     pause(1) % Pause for 1 seconds (unit)
end
delete(t)

CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201]; % Wavelength range (380-780 nm) with 2 nm interval

% Set random channel combinations to turn on
numsample = 30; % Set the number of the samples to be measured
SpectrumInputSet = [0,150,200,252]; % Set the current set to be randomized (0=off)

% Generate the random combinations of the LED channels
for i = 1:numsample
    SpectrumInput(i,:) = datasample(SpectrumInputSet,16); % Set current randomly per each channel from the 'current_set_combo'
end

% Temp for measurements
SpectrumInput(:,1:5) = 0;
SpectrumInput(:,9) = 0;

% Vpixx control
for t = 1:numsample    
   
    % Set SubColors
    for k=1:16
        Datapixx('SetPropixxHSLedCurrent', 0, k-1, SpectrumInput(t,k)); % Primary 1
        Datapixx('SetPropixxHSLedCurrent', 1, k-1, SpectrumInput(t,k)); % Primary 2
        Datapixx('SetPropixxHSLedCurrent', 2, k-1, SpectrumInput(t,k)); % Primary 3
    end
   
    % PR670 measure
    fw_rand(:,t)=MeasSpd(S,5,'all');
    disp(append('Measurement ',num2str(t),' complete!'))
end

% Save the spectrums (.mat file)
filename = append(Date,'SpectrumInput','.mat');
save(filename,'SpectrumInput')

filename = append(Date,'Spectrum_rand','.mat');
save(filename,'fw_rand')

%% Plot Power spectrum
figure(1); subplot(2,2,1); hold on;
plot(SToWls(S),fw,'b-');
xlabel('Wavelength(nm)')
ylabel('Spectral irradiance')
xlim([380 780]);

% CIE (x,y) chromaticity
figure(1);subplot(2,2,2); hold on;
plot(xyY(1,:),xyY(2,:),'r.'); % Measurement point
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
legend('Test');
