function [targetSpdMeasured] = MeasureDesiredTargetPrimaries(targetPrimaries,subprimaryNInputLevels,subPrimaryCalstructData,options)
% Measure the desired target primaries to use them for computing the contrast image.
%
% Syntax:
%    [targetSpdMeasured] = MeasureDesiredTargetPrimaries(targetPrimaries,subprimaryNInputLevels,subPrimaryCalstructData,options)
%
% Description:
%    This measures the target primaries that we actually get and then use
%    the measured rather than the nominal primaries to compute the image.
%
% Inputs:
%    targetPrimaries -            Target primaries you wish to reproduce.
%    subprimaryNInputLevels -     Device max input levels for the
%                                 subprimaries.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with.
%
% Outputs:
%    targetSpdMeasured -          Measured SPDs of the desired isolating
%                                 target primaries.
%
% Optional key/value pairs:
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either 'Normal' (true) or
%                                 'Steady-on' (false).
%
% History:
%    10/06/21  smo  Started on it

%% Set parameters.
arguments
    targetPrimaries
    subprimaryNInputLevels
    subPrimaryCalstructData
    options.projectroMode (1,1) = true
end

%% Initialize (connect to both display and measurement device).
%
% This part prepares for the measurements which includes connecting PC to the
% projector, setting the projector display state (normal/steady-on modes),
% and talking to the measurement device to ready to run (PR670).
%
% Add VPixx toolbox 'Datapixx' to the path.
toolboxDirectory = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx'; % This is where Linux box store the file.
addpath(genpath(toolboxDirectory));

% Connect to the Vpixx projector.
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

% Set the projector mode.
if (options.projectroMode)
    commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Normal mode (Default)
    unix(commandNormal)
    disp('Projector is set as Normal mode');
else
    commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit'; % Steady-on mode
    unix(commandSteadyOn)
    disp('Projector is set as Steady-on mode');
end

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

% Set the working range of subprimary channels (chanel 8 is not working at the moment).
logicalToPhysical = [0:7 9:15];

% Measurement range.
S = subPrimaryCalstructData.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalstructData.get('nDevices');

% Target Spd.
% This will be compared to the measurement result at the end.
for pp = 1:nPrimaries
    targetSpd(:,pp) = PrimaryToSpd(subPrimaryCalstructData,targetPrimaries(:,pp));
end

%% Set primary settings and measure it.
% Display plain screen on DLP using PTB.
PsychDefaultSetup(2); % PTB pre-setup
screens = Screen('Screens');
screenNumber = max(screens);
white = WhiteIndex(screenNumber);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
%backgroundColor = [255 255 255];
%[window, backgroundRect] = Screen('OpenWindow', screenNumber, backgroundColor, windowRect);

% Set projector settings and measure.
otherPrimarySettings = 0;
for pp = 1:nPrimaries
    otherPrimaries = setdiff(1:nPrimaries,pp);
    % Set target settings.
    targetSettings = PrimaryToSettings(subPrimaryCalstructData,targetPrimaries((:,pp)) * subprimaryNInputLevels;
    for ss = logicalToPhysical(1:nSubprimaries)
        % Set projector current levels.
        Datapixx('SetPropixxHSLedCurrent', pp-1, ss, round(targetSettings(ss))); % Target primary
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(1)-1, ss, otherPrimarySettings); % Other Primary 1
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(2)-1, ss, otherPrimarySettings); % Other Primary 2
    end
    % Measurement.
    targetSpdMeasured(:,pp) = MeasSpd(S,5,'all');
end

% Close PTB screen.
sca;

% Plot the measurement result.
figure;
for pp = 1:nPrimaries;
    subplot(1,pp);
    plot(SToWls(S),targetSpd(:,pp),'k-','LineWidth',3); % Target Spd.
    plot(SToWls(S),targetSpdMeasured(:,pp),'r--','LineWidth',3); % Measured Spd.
    xlabel('Wavelength (nm)');
    ylabel('Spectral power');
    legend('Target','Measurement');
    title(sprintf('Primary %d',pp));
end

end