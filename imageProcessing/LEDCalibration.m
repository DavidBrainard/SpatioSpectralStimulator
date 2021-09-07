%% Calibration
% Data will be saved at the folder (PsychCallDemoData)

%% Initialize 
clear all; close all; clc;
tbUse('BrainardLabBase')
tbUseProject('SpatioSpectralStimulator'); % Hook to the Dropbox for saving the calibration file

addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx')); % add VPixx toolbox ('Datapixx') to path (incl. all subfolders)

% Connect to PR-670
command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
command_check = 'ls -l /dev/ttyACM0';
unix(command_check)
command_check = 'ls -l /dev/ttyACM1';
unix(command_check)

% Connect to the Vpixx projector
command_primaries = 'vputil rw 0x1c8 0x7 -q quit'; % Set to 'All Primaries on'
unix(command_primaries)

isReady = Datapixx('open');
isReady = Datapixx('IsReady'); 

%% Vpixx projector color control (useful for warming up the device)
current = 252;
for j=1:3 % PrimaryColor
   for i=1:16 % SubColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current); 
   end
end

for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 
currents

%% Set primries for calibration
% ***************Control by each primary color channel***************
% First, set all currents to baseline
isReady = Datapixx('open'); 
isReady = Datapixx('IsReady');

current = 0;

for j=1:3 % PrimaryColor
   for i=1:16 % SubColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current); 
   end
end

% Set one peak per each primary for calibration
current = 252;

subcolor_R = 13; % 0-15 / 8 is not working
subcolor_G = 7; % 0-15
subcolor_B = 2; % 0-15

Datapixx('SetPropixxHSLedCurrent',0,subcolor_R,current) % R
Datapixx('SetPropixxHSLedCurrent',1,subcolor_G,current) % G
Datapixx('SetPropixxHSLedCurrent',2,subcolor_B,current) % B

for j=1:3; 
    for i=1:16; 
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

currents

%% Execution part
% [0.6] meters between DLP and PR670
OOC_calibrateMonitor

%% Analyze
OOC_analyzeCal

%% How to use 
% Hi Semin,
% I have added a case for the SACC display in OOC_CalibrateMonitor.
% 
% 0. Please pull the BrainardLabToolbox repository to get the latest version.
% 00. Make sure the PR650 is setup, placed 1 m from the center of the display and that you can talk to it.
% 
% The settings for SACC are in method generateConfigurationForSACC() which starts in line 576. 
% 
% You may need to change
% displaySettings.screenToCalibrate, depending on whether the SACC display is set up as the main screen or the secondary screen. Right now, the code sets it to 1, which assumes the SACC is the main screen.
% 
% Also you may want to change some CalibratorOptions, to speed up the first calibration:
% For example, the ’nMeas’, which is the number of samples along the gamma curve, which are now set to 21, or'nAverage’, which is the number of repeated measurements which are now set to 2.
% 
% 
% 1. Once you have done some (or none) of the above changes, run OOC_CalibrateMonitor.
% You will be prompted with:
% 
% Available calibration configurations 
% MetropsisCalibration
% VirtualWorldCalibration
% NaturaImageThresholds
% ColorMaterialCalibration
% ViewSonicProbe
% BOLDscreen
% SamsungOLEDpanel
% Left_SONY_PVM2541A
% Right_SONY_PVM2541A
% InVivoSensaVue_FlatPanel
% AppleThunderboltDisplay
% HTCVive
% SACC
% Select a calibration config [NaturaImageThresholds]: 
% 
% Enter SACC
% 
% 2. Then the program will ask for you name: Enter it.
% 
% 3. Following this, the program will ask you for an email to send the notification when the calibration is done. Right now it is set to your sas upenn email. You can enter a different one, if you wish so. Or hit return to accept the default email.
% 
% 4. Then you will be prompted with
% Available radiometer types:
% [  1]. PR650dev
% [  2]. PR670dev
% [  3]. SpectroCALdev
% Select a radiometer type (1-3) [1]: 
% 
% Enter the one you will be using, 1 if it is the PR650.
% 
% 5. Finally, you will be prompted with:
% Available calibrator types:
%         [ 1]. MGL-bases
%         [ 2]. PsychImaging-based (8-bit).
%         Select a calibrator type (1-2) [1]
% 
% Enter 2, to select the PTB PsychImaging based engine.
% 
% Then follow the instructions to start the calibration. At some point the program will tell you that you have 3 seconds to leave the room. Assuming all lights are off and everything else set up, hit return and leave the room. You will get an email when the calibration is done, which will probably be 1 hour.
% 
% The first time you can stay in the room to make sure the program has started and that there is no crash.
% 
% You can view the results of the calibration by running OOC_analyzeCal.
% 
% Good luck, and let me know if there is a problem.
% 
% Best,
% Nicolas

%% Just for my record

% % Configuration function for the SACC display (LED/DLP optical system)
% function [displaySettings, calibratorOptions] = generateConfigurationForSACC()
%     % Specify where to send the 'Calibration Done' notification email
%     emailAddressForNotification = 'seminoh@sas.upenn.edu';
%     
%     % Specify the @Calibrator's initialization params. 
%     % Users should tailor these according to their hardware specs. 
%     % These can be set once only, at the time the @Calibrator object is instantiated.
%     displaySettings = { ...
%         'screenToCalibrate',        2, ...                          % which display to calibrate. main screen = 1, second display = 2
%         'desiredScreenSizePixel',   [1920 1080], ...                % pixels along the width and height of the display to be calibrated
%         'desiredRefreshRate',       120, ...                         % refresh rate in Hz
%         'displayPrimariesNum',      3, ...                          % for regular displays this is always 3 (RGB) 
%         'displayDeviceType',        'monitor', ...                  % this should always be set to 'monitor' for now
%         'displayDeviceName',        'SACC', ...                     % a name for the display been calibrated
%         'calibrationFile',          'SACC', ...                     % name of calibration file to be generated
%         'comment',                  'The SACC LED/DLP optical system' ...          % some comment, could be anything
%         };
%     
%     % Specify the @Calibrator's optional params using a CalibratorOptions object
%     % To see what options are available type: doc CalibratorOptions
%     % Users should tailor these according to their experimental needs.
%     calibratorOptions = CalibratorOptions( ...
%         'verbosity',                        2, ...
%         'whoIsDoingTheCalibration',         input('Enter your name: ','s'), ...
%         'emailAddressForDoneNotification',  GetWithDefault('Enter email address for done notification',  emailAddressForNotification), ...
%         'blankOtherScreen',                 0, ...                          % whether to blank other displays attached to the host computer (1=yes, 0 = no), ...
%         'whichBlankScreen',                 1, ...                          % screen number of the display to be blanked  (main screen = 1, second display = 2)
%         'blankSettings',                    [0.0 0.0 0.0], ...              % color of the whichBlankScreen 
%         'bgColor',                          [0.3962 0.3787 0.4039], ...     % color of the background  
%         'fgColor',                          [0.3962 0.3787 0.4039], ...     % color of the foreground
%         'meterDistance',                    1.0, ...                        % distance between radiometer and screen in meters
%         'leaveRoomTime',                    3, ...                          % seconds allowed to leave room
%         'nAverage',                         2, ...                          % number of repeated measurements for averaging
%         'nMeas',                            21, ...                         % samples along gamma curve
%         'boxSize',                          600, ...                        % size of calibration stimulus in pixels (it was 150 / Semin)
%         'boxOffsetX',                       0, ...                          % x-offset from center of screen (neg: leftwards, pos:rightwards)         
%         'boxOffsetY',                       0 ...                           % y-offset from center of screen (neg: upwards, pos: downwards)                      
%     );
% 
% end


