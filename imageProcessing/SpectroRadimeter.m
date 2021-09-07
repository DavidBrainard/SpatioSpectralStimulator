%% VPixx - setting the ground state of the projector

% Initialize
clear all; close all; clc;
tbUse('BrainardLabBase')

%% Measurement using PR-670 Spectroradiometer
% In terminal, command to authorize to read and write 'sudo chmod a+rw /dev/ttyACM0 (the last character is number zero)'
% Also command to check the status 'ls -l /dev/ttyACM0'

%% Initialize (authorization)
% Mostly PR670 connects to 'ttyACM0', but sometimes to 'ttyACM1', so just
% do both to make it sure to be connected

command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)

command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)

command_check = 'ls -l /dev/ttyACM0';

unix(command_check)

%% Measurement 
CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670

S=[380 2 201]; % Wavelength range (380-780 nm) / S=[380 2 201]
fw=MeasSpd(S,5,'all'); % Measurement (it takes a while)
plot(SToWls(S),fw,'b-');

%% Save the spectrum

filename = append('0729Black30(ch6-15)','.mat'); 
save(filename,'fw')

%% Further calculations
load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw'*T_XYZ; % XYZ calculation
xyY = XYZToxyY(XYZ'); % xyY
colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

XYZw = [100 100 100]'; % Lab white
Lab = XYZToLab(XYZ',XYZw); % Lab calculation

% Power spectrum
figure(1); subplot(2,2,1); hold on;
plot(SToWls(S),fw,'b-');
xlabel('Wavelength(nm)')
ylabel('Spectral irradiance')
xlim([380 780]);

% CIE (x,y) chromaticity
figure(1);subplot(2,2,2); hold on;
plot(xyY(1,:),xyY(2,:),'r.'); % Measurement point
plot(colorgamut(1,:),colorgamut% Configuration function for the SACC display (LED/DLP optical system)
function [displaySettings, calibratorOptions] = generateConfigurationForSACC()
    % Specify where to send the 'Calibration Done' notification email
    emailAddressForNotification = 'seminoh@sas.upenn.edu';
    
    % Specify the @Calibrator's initialization params. 
    % Users should tailor these according to their hardware specs. 
    % These can be set once only, at the time the @Calibrator object is instantiated.
    displaySettings = { ...
        'screenToCalibrate',        2, ...                          % which display to calibrate. main screen = 1, second display = 2
        'desiredScreenSizePixel',   [1920 1080], ...                % pixels along the width and height of the display to be calibrated
        'desiredRefreshRate',       120, ...                         % refresh rate in Hz
        'displayPrimariesNum',      3, ...                          % for regular displays this is always 3 (RGB) 
        'displayDeviceType',        'monitor', ...                  % this should always be set to 'monitor' for now
        'displayDeviceName',        'SACC', ...                     % a name for the display been calibrated
        'calibrationFile',          'SACC', ...                     % name of calibration file to be generated
        'comment',                  'The SACC LED/DLP optical system' ...          % some comment, could be anything
        };
    
    % Specify the @Calibrator's optional params using a CalibratorOptions object
    % To see what options are available type: doc CalibratorOptions
    % Users should tailor these according to their experimental needs.
    calibratorOptions = CalibratorOptions( ...
        'verbosity',                        2, ...
        'whoIsDoingTheCalibration',         input('Enter your name: ','s'), ...
        'emailAddressForDoneNotification',  GetWithDefault('Enter email address for done notification',  emailAddressForNotification), ...
        'blankOtherScreen',                 0, ...                          % whether to blank other displays attached to the host computer (1=yes, 0 = no), ...
        'whichBlankScreen',                 1, ...                          % screen number of the display to be blanked  (main screen = 1, second display = 2)
        'blankSettings',                    [0.0 0.0 0.0], ...              % color of the whichBlankScreen 
        'bgColor',                          [0.3962 0.3787 0.4039], ...     % color of the background  
        'fgColor',                          [0.3962 0.3787 0.4039], ...     % color of the foreground
        'meterDistance',                    1.0, ...                        % distance between radiometer and screen in meters
        'leaveRoomTime',                    3, ...                          % seconds allowed to leave room
        'nAverage',                         2, ...                          % number of repeated measurements for averaging
        'nMeas',                            21, ...                         % samples along gamma curve
        'boxSize',                          600, ...                        % size of calibration stimulus in pixels (it was 150 / Semin)
        'boxOffsetX',                       0, ...                          % x-offset from center of screen (neg: leftwards, pos:rightwards)         
        'boxOffsetY',                       0 ...                           % y-offset from center of screen (neg: upwards, pos: downwards)                      
    );

end(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
legend('Test');

% CIELAB 
figure(1);subplot(2,2,3); hold on;
plot(Lab(1,:),Lab(2,:),'r.'); % Measurement point
xlabel('CIELAB a*')
ylabel('CIELAB b*')
xlim([-100 100]);
ylim([-100 100]);
legend('Test');


%% Cone signals calculation
% load T_cones_sp % Smith-Pokorny Cone spectral sensitivity function
% T_Cones = T_cones_sp';
% Cones = fw'.*T_Cones; % Cone signals calculation

