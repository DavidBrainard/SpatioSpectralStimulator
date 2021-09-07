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
% The password should be typed at the very first command

command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
command_check = 'ls -l /dev/ttyACM0';
unix(command_check)
command_check = 'ls -l /dev/ttyACM1';
unix(command_check)

% Start 2 (Vpixx projector)
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

%% Vpixx projector color control (useful for warming up the device)

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

%% Single measurement (TEST)
CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670

S=[380 2 201]; % Wavelength range (380-780 nm) / S=[380 2 201]
fw=MeasSpd(S,5,'all'); % Measurement (it takes a while)
plot(SToWls(S),fw,'b-');

%% ******************************************************ROUTINE STARTS FROM HERE****************************************************************
%% SINGLE CHANNEL MEASUREMENT
% Measurement time = 18 min / 68 test colors (48 singles + 20 randoms)
clear start; clear fw_single; clear fw_white; clear fw_blk; clear fw_rand; clear SpectrumInput;

CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201]; % Wavelength range (380-780 nm) with 2 nm interval
Date = 'Dataset3-3)';
ADD = '(WithinPrimary)_black10';

% Black and white levels
current_white = 252; 
current_blk = 10; % Set black level (10 is measurable)
current_base = 10; % Set the black level as a base
currentset = [252];
datatype = 1; % 1 = within primary / 2 = across primary
numsample = 20; % Number of random spectra to generate

% Input data needs to be filled till this part
% *************************************************************************

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

% White measurement
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

for j=1:3 % PrimaryColor
   for i=1:16 % subColor             
       Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current_white); 
   end
end

fw_white = MeasSpd(S,5,'all'); 
disp(append('White ','measurement',' complete!'))

% Black measurement
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

for j=1:3 % PrimaryColor
   for i=1:16 % subColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current_blk); 
   end
end

fw_blk = MeasSpd(S,5,'all');
disp(append('Black ','measurement',' complete!'))

% Single peak spectrum measurements
if datatype == 1

    % 1) WITHIN PRIMARY (All primaries have a same value)
    for i=1:16 % Num of 15 LED channels (i=1-16 / Vpixx input=0-15)

        % Set all channel currents to 0 as baseline
        isReady = Datapixx('open');
        isReady = Datapixx('IsReady');

        for j=1:3 % PrimaryColor
            for k=1:16 % subColor
             Datapixx('SetPropixxHSLedCurrent', j-1, k-1, current_base);
            end
        end

          % currentset = [10:10:240,252]; % p.s. These all sets took 1hr 45min (a total of 368 measurements)
          currentset = [252];

          for f = 1:length(currentset);
            current = currentset(f);

            % Measure each LED channel
            isReady = Datapixx('open');
            isReady = Datapixx('IsReady');

            Datapixx('SetPropixxHSLedCurrent', 0, i-1, current); % Primary 1
            Datapixx('SetPropixxHSLedCurrent', 1, i-1, current); % Primary 2
            Datapixx('SetPropixxHSLedCurrent', 2, i-1, current); % Primary 3

            % PR670 measurement
            fw_single(:,f+length(currentset)*(i-1))=MeasSpd(S,5,'all');
            disp(append('LED Channel ',num2str(i),' current ',num2str(f),' complete!'))
          end
    end
    
elseif datatype == 2
    
    % 2) ACROSS PRIMARY (Primary has a different value)
    for p=1:3 % Primary

          for s = 1:16 % Subprimary

                % Set all channel currents as baseline
                isReady = Datapixx('open');
                isReady = Datapixx('IsReady');

                for j=1:3 % PrimaryColor
                    for k=1:16 % subColor
                     Datapixx('SetPropixxHSLedCurrent', j-1, k-1, current_base);
                    end
                end

                for c = 1:length(currentset);
                    current = currentset(c);

                    % Measure each LED channel
                    isReady = Datapixx('open');
                    isReady = Datapixx('IsReady');
                    Datapixx('SetPropixxHSLedCurrent', p-1, s-1, current);

                    % PR670 measurement
                    fw_single(:,c + length(currentset)*(s-1) + 16*length(currentset)*(p-1))=MeasSpd(S,5,'all');
                    disp(append('Primary ',num2str(p),' Subchannel ',num2str(s),' current ',num2str(c),' complete!'))
                end
          end
      
    end
end

% Save the spectrums (.mat file)
filename = append(Date,'Spectrum_single',ADD,'.mat');            
save(filename,'fw_single') 

filename2 = append(Date,'Spectrum_white',ADD,'.mat');            
save(filename2,'fw_white')

filename3 = append(Date,'Spectrum_black',ADD,'.mat');        
save(filename3,'fw_blk') 

% ******************************************************ROUTINE STARTS FROM HERE****************************************************************
% RANDOM COMBO OF SPECTRUM
clear start; clear SpectrumInput; clear fw_rand;

% Date = '0805';
% ADD = '(AcrossPrimary)_black0';

% Time delay before start (go out the room before it's done)
% timesec = 15;
% t = timer('TimerFcn','stat=false; disp(''Time OFF to leave!'')','StartDelay',timesec);
% start(t)
% stat=true;
% while(stat==true)
%      disp('.')
%      pause(1) % Pause for 1 seconds (unit)
% end
% delete(t)

CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201]; % Wavelength range (380-780 nm) with 2 nm interval

% % Generate random combinations
% numsample = 20; % Set the number of the samples to be measured
%     
% if datatype == 1
%     
%     % 1) WITHIN PRIMARY
% %     SpectrumInput = ones(numsample,16).*current_base; % Base structure
%     SpectrumInputSet = [current_base,252]; % Set the current set to be randomized (current_base=off)
% 
%     for i = 1:numsample
%         SpectrumInput(i,:) = datasample(SpectrumInputSet,16); % Set current randomly per each channel from the 'currentset_combo'
%     end
% 
%     SpectrumInput(:,9) = 0; % Turn off for specific LED Channels (Ch 9 not working)
%     
% elseif datatype == 2
%     
%     % 2) ACROSS PRIMARY
%     SpectrumInput = ones(numsample*3,16).*current_base; % Base structure
%     
%     SubColorSet = [1:8,10:16]; % Ch 1-16 without Ch 9 which is not working
% 
%     for i = 1:numsample*3
%         idx_RandSubColor(i) = datasample(SubColorSet,1); % Set SubColor channel (rand)
%         idx_SpectrumInput(i) = datasample(CurrentSet,1); % Set current value (rand)
%     end
% 
%     for i = 1:numsample*3
%         SpectrumInput(i,idx_RandSubColor(i))=idx_SpectrumInput(i);
%     end
% 
%     SpectrumInput_R = SpectrumInput(1:numsample,:); % Primary 1
%     SpectrumInput_G = SpectrumInput(1+numsample:numsample*2,:); % Primary 2
% %     SpectrumInput_B = SpectrumInput(1+numsample*2:numsample*3,:); % Primary 3
%     
%     SpectrumInput_B = ones(numsample,16).* current_base; % Primary 3 turnoff
%     
%     SpectrumInput = {SpectrumInput_R SpectrumInput_G SpectrumInput_B}; % SpectrumInput {R G B}
% end


% Load the pre-saved SpectrumInput set
load('Dataset3)SpectrumInput(WithinPrimary)_black10'); % 'SpectrumInput'

% Correct current base level on loaded SpectrumInput
% SpectrumInput_base = ones(numsample,16).*current_base;
% 
% if datatype == 1
%     SpectrumInput = SpectrumInput_base + max(SpectrumInput - SpectrumInput_base,0);
% elseif datatype == 2
%     SpectrumInput_R = SpectrumInput_base + max(cell2mat(SpectrumInput(1)) - SpectrumInput_base,0); % Primary 1
%     SpectrumInput_G = SpectrumInput_base + max(cell2mat(SpectrumInput(2)) - SpectrumInput_base,0); % Primary 2
%     SpectrumInput_B = SpectrumInput_base + max(cell2mat(SpectrumInput(3)) - SpectrumInput_base,0); % Primary 3
%     SpectrumInput = {SpectrumInput_R SpectrumInput_G SpectrumInput_B};
% end

% VPixx project current control
if datatype == 1

    % 1) WITHIN PRIMARY (All Primaries have the same value)current = 252; % 0-252

    for t = 1:numsample    

        % Set SubColors
        for k=1:16
            Datapixx('SetPropixxHSLedCurrent', 0, k-1, SpectrumInput(t,k)); % Primary 1
            Datapixx('SetPropixxHSLedCurrent', 1, k-1, SpectrumInput(t,k)); % Primary 2
            Datapixx('SetPropixxHSLedCurrent', 2, k-1, SpectrumInput(t,k)); % Primary 3
        end

        % PR670 measure
        fw_rand(:,t)=MeasSpd(S,5,'all');
        disp(append('Rand Measurement ',num2str(t),' complete!'))
    end

elseif datatype == 2

    % 2) ACROSS PRIMARY (Only one primary is turned on)
    for t = 1:numsample    

          % Set SubColors with the different value
        for k=1:16
            Datapixx('SetPropixxHSLedCurrent', 0, k-1, SpectrumInput_R(t,k)); % Primary 1
            Datapixx('SetPropixxHSLedCurrent', 1, k-1, SpectrumInput_G(t,k)); % Primary 2
            Datapixx('SetPropixxHSLedCurrent', 2, k-1, SpectrumInput_B(t,k)); % Primary 3
        end

           % PR670 measure
            fw_rand(:,t)=MeasSpd(S,5,'all');
            disp(append('Testcolor ',num2str(t),' measurement complete!'))
    end
end 

% Save the spectrums (.mat file)
filename = append(Date,'SpectrumInput',ADD,'.mat');
save(filename,'SpectrumInput')

filename = append(Date,'Spectrum_rand',ADD,'.mat');
save(filename,'fw_rand')


%% ******************************************************ROUTINE STARTS FROM HERE****************************************************************
%% STABILITY TEST OVER TIME
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
