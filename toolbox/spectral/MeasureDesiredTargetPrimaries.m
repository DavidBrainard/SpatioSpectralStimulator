function [targetSpdMeasured] = MeasureDesiredTargetPrimaries(targetPrimaries,subPrimaryCalStructData, ...
    targetPrimaryNum,options)
% Measure the desired target primaries to use them for computing the contrast image.
%
% Syntax:
%    [targetSpdMeasuredRaw,targetSpdMeasuredNorm,measuredToTarget] = MeasureDesiredTargetPrimaries(targetPrimaries,subPrimaryCalStructData, ...
%                                                                    targetPrimaryNum,options)
%
% Description:
%    This measures the target primaries that we actually get and then use
%    the measured rather than the nominal primaries to compute the image.
%
% Inputs:
%    targetPrimaries -            Target primaries you wish to reproduce.
%                                 These have not yet been gamma corrected.
%    subPrimaryCalStructData -    The calibration object that describes the
%                                 device we're working with.
%    targetPrimaryNum -           Target primary number. This is required
%                                 when setting up the projector current input
%                                 levels.
%
% Outputs:
%    targetSpdMeasured -          Measured SPDs of the desired isolating
%                                 target primaries.
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
%    'verbose'           -        Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    10/06/21  smo                Started on it
%    10/07/21  smo                Makes it working for a single primary and
%                                 added a feature to skip the measurement part.
%    10/19/21  smo                Added to save out the normalized spd results.
%    10/20/21  smo                Discard the normalized spd outputs.

%% Set parameters.
arguments
    targetPrimaries
    subPrimaryCalStructData
    targetPrimaryNum
    options.projectorMode (1,1) = true
    options.measurementOption (1,1) = true
    options.verbose (1,1) = true
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

% Get number of discrete input levels for the device.  
% So, control values go from 0 to (subprimaryNInputLevels-1).
subprimaryNInputLevels = size(subPrimaryCalStructData.get('gammaInput'),1);

% Measurement range.
S = subPrimaryCalStructData.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalStructData.get('nDevices');
if (targetPrimaryNum > 3) && (targetPrimaryNum < 1)
    error('Target primary number should be set within [1, 2, 3]')
end

% Target Spd.
% This will be compared to the measurement result at the end.
targetSpd = PrimaryToSpd(subPrimaryCalStructData,targetPrimaries);

%% Set primary settings and measure it.
if (options.measurementOption)
    % Display plain screen on DLP using PTB.
    PsychDefaultSetup(2); % PTB pre-setup
    screens = Screen('Screens');
    screenNumber = max(screens);
    white = WhiteIndex(screenNumber);
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
    
    % Set projector input settings.  These are gamma corrected but still
    % live as real numbers on the interval [0,1], following the convention
    % we use for other types of display.
    targetSubprimarySettings = PrimaryToSettings(subPrimaryCalStructData,targetPrimaries);
    otherPrimaries = setdiff(1:nPrimaries,targetPrimaryNum);
    otherPrimarySubprimarySettings = 0;
    
    % Set projector current levels as the above settings.
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', targetPrimaryNum-1, logicalToPhysical(ss), round(targetSubprimarySettings(ss)*(subprimaryNInputLevels-1)));   % Target primary
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(1)-1, logicalToPhysical(ss), round(otherPrimarySubprimarySettings*(subprimaryNInputLevels-1))); % Other Primary 1
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(2)-1, logicalToPhysical(ss), round(otherPrimarySubprimarySettings*(subprimaryNInputLevels-1))); % Other Primary 2
    end
    
    % Measurement.
    targetSpdMeasured = MeasSpd(S,5,'all');
    if (options.verbose)
        fprintf('Measurement complete! - Primary %d\n',targetPrimaryNum);
    end
else
    targetSpdMeasured = targetSpd; % Just put the same spd as desired for here.
    if (options.verbose)
        fprintf('Measurement has been skipped!\n');
    end
end

% Close PTB screen.
sca;

% Plot the results.
if (options.verbose)
    figure; hold on;
    plot(SToWls(S),targetSpd,'k-','LineWidth',1); % Target Spd.
    plot(SToWls(S),targetSpdMeasured,'r--','LineWidth',1); % Measured raw spd.
    xlabel('Wavelength (nm)');
    ylabel('Spectral power');
    legend('Target','Measurement_Raw','Meausurement_Norm');
    title(sprintf('Primary %d',targetPrimaryNum));
end
end