function [measuredLMSContrast] = MeasureLMSContrastGaborPatch(targetLMSContrast,bgPrimaries,subprimaryNInputLevels,subPrimaryCalstructData,targetPrimaryNum,options)
% Measure the LMS contrasts at some points on the gabor patch image.
%
% Syntax:
%
%
% Description:
%    This measures and calculates the LMS contrasts of the modulated gabor
%    patch image, so we can check if the modulation was done properly.
%
% Inputs:
%    targetLMSContrast -          To be filled
%    subprimaryNInputLevels -     Device max input levels for the
%                                 subprimaries.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with.
%    targetPrimaryNum -           Target primary number. This is required
%                                 when setting up the projector current input
%                                 levels.
%
% Outputs:
%    measuredLMSContrast -        Measured LMS contrast of the modulated
%                                 gabor patch image.
%
% Optional key/value pairs:
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either 'Normal' (true) or
%                                 'Steady-on' (false).
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%
% History:
%    10/08/21  smo                Started on it

%% Set parameters.
arguments
    subprimaryNInputLevels
    subPrimaryCalstructData
    targetPrimaryNum
    options.projectorMode (1,1) = true
    options.measurementOption (1,1) = true
end

%% Initialize (connect to both display and measurement device).
%
% This part prepares for the measurements which includes connecting PC to the
% projector, setting the projector display state (normal/steady-on modes),
% and talking to the measurement device to ready to run (PR670).
%
% Measurement pre-preparations.
if (options.measurementOption)
    % Connect to the Vpixx projector.
    isReady = Datapixx('open');
    isReady = Datapixx('IsReady');
    
    % Add VPixx toolbox 'Datapixx' to the path.
    toolboxDirectory = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx'; % This is where Linux box store the file.
    addpath(genpath(toolboxDirectory));
    
    % Give read-and-write permission to the spectroradiometer.
    % This part is exclusivley when using the Linux box, and
    % the PR670 is connected to the device either 'ttyACM0' or 'ttyACM1',
    % so give permission to both here.
    commandConnectToPR670_1 = 'sudo chmod a+rw /dev/ttyACM0';
    unix(commandConnectToPR670_1)
    commandConnectToPR670_2 = 'sudo chmod a+rw /dev/ttyACM1';
    unix(commandConnectToPR670_2)
    
    % Connect to spectroradiometer PR670.
    CMCheckInit(5); % Number 5 is allocated to PR670
    
    % Set the projector mode.
    if (options.projectorMode)
        commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Normal mode (Default)
        unix(commandNormal)
        disp('Projector is set as Normal mode');
    else
        commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit'; % Steady-on mode
        unix(commandSteadyOn)
        disp('Projector is set as Steady-on mode');
    end
end

% Set the working range of subprimary channels (chanel 8 is not working at the moment).
logicalToPhysical = [0:7 9:15];

% Measurement range.
S = subPrimaryCalstructData.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalstructData.get('nDevices');
if (targetPrimaryNum > 3) && (targetPrimaryNum < 1)
    error('Target primary number should be set within [1, 2, 3]')
end

%% Set the desired LMS contrast point and measure it. This part should be updated later.
if (options.measurementOption)
    % Display plain screen on DLP using PTB.
    PsychDefaultSetup(2); % PTB pre-setup
    screens = Screen('Screens');
    screenNumber = max(screens);
    white = WhiteIndex(screenNumber);
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
    
    % Set projector input settings.
    targetSettings = PrimaryToSettings(subPrimaryCalstructData,targetPrimaries) * subprimaryNInputLevels;
    otherPrimaries = setdiff(1:nPrimaries,targetPrimaryNum);
    otherPrimarySettings = 0;
    
    % Set projector current levels as the above settings.
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', targetPrimaryNum-1, logicalToPhysical(ss), round(targetSettings(ss))); % Target primary
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(1)-1, logicalToPhysical(ss), otherPrimarySettings); % Other Primary 1
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(2)-1, logicalToPhysical(ss), otherPrimarySettings); % Other Primary 2
    end
    
    % Measurement.
    targetSpdMeasured = MeasSpd(S,5,'all');
    disp(sprintf('Measurement complete! - Primary %d',targetPrimaryNum));
else
    targetSpdMeasured = targetSpd; % Just put the same spd as desired for here.
    disp(sprintf('Measurement has been skipped!'));
end

% Close PTB screen.
sca;



%% Measure some points on the Gabor patch for calculating LMS contrast.


end
