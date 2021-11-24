% SACC_ScreenStabilityCheck
%
% This tests screen stability for SACC project. We want to know how 
% 
% History:
%    11/24/2021 smo   Pulled out the part from the old code and cleaned
%                     tried to follow our convention for variable names etc.

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