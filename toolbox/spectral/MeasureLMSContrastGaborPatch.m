function [targetLMSContrastMeasured] = MeasureLMSContrastGaborPatch(contrastImage,isolatingPrimaries,projectorCalObj,bgSpd,...
    subprimaryNInputLevels,subPrimaryCalstructData,options)
% Measure the LMS contrasts at some points on the gabor patch image.
%
% Syntax: [targetLMSContrastMeasured] = MeasureLMSContrastGaborPatch(quantizedContrastImage,isolatingPrimaries,projectorCalObj,bgSpd,...
%    subprimaryNInputLevels,subPrimaryCalstructData)
%
%
% Description:
%    This measures and calculates the LMS contrasts of the modulated gabor
%    patch image, so we can check if the modulation was done properly.
%
% Inputs:
%    contrastImage -              Modulated contrast gabor patch image.
%    isolatingPrimaries -         Isolating primaries which
%                                 reproduce the desired LMS contrasts.
%    projectorCalObj -            The calibration object that contains the
%                                 the data for the projector.
%    bgSpd                        Background Spd of the gabor patch.
%    subprimaryNInputLevels -     Device max input levels for the
%                                 subprimaries.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with.
%
% Outputs:
%    targetLMSContrastMeasured -  Measurement results of the some points on
%                                 the contrast gabor patch image.
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
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    10/08/21  smo                Started on it
%    10/14/21  smo                Made a working draft.

%% Set parameters.
arguments
    contrastImage
    isolatingPrimaries
    projectorCalObj
    bgSpd
    subprimaryNInputLevels
    subPrimaryCalstructData
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

% Measurement range.
S = subPrimaryCalstructData.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalstructData.get('nDevices');

%% Pick some points on the Gabor patch and calculate target settings for measurements.
bgLMS = T_cones * sum(bgSpd,2);
imageN = size(contrastImage,1);
centerN = imageN/2;

% Set measurement points on the Gabor patch. Points are decided on the
% X-pixel 2-D plane.
measureXPixelPoints = round([0.3:0.1:0.7] * imageN);
nTestPoints = size(measureXPixelPoints,2);

% Calculate LMS contrast.
SetGammaMethod(projectorCalObj,0);
measureDesiredLMSContrast = squeeze(contrastImage(centerN,measureXPixelPoints,:))';
measureDesiredLMS = ContrastToExcitation(measureDesiredLMSContrast,bgLMS);
measureDesiredPrimaries = SensorToPrimary(projectorCalObj,measureDesiredLMS);

% Get target settings for the projector.
for tt = 1:nTestPoints
    % Set target primaries.
    measureTargetPrimaries(:,:,tt) = isolatingPrimaries .* measureDesiredPrimaries(:,tt)';
    for pp = 1:nPrimaries
        % Set target settings per each test point.
        measureTargetSettings(:,pp,tt) = PrimaryToSettings(subPrimaryCalstructData{pp},measureTargetPrimaries(:,pp,tt));
    end
end

%% Measure it.
% Set projector setting for each test point on the gabor patch and measure it.
if (options.measurementOption)
    % Display plain screen on DLP using PTB.
    PsychDefaultSetup(2); % PTB pre-setup
    screens = Screen('Screens');
    screenNumber = max(screens);
    white = WhiteIndex(screenNumber);
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
    
    % Set projector current levels as the above settings.
    for tt = 1:nTestPoints
        % Set projector current levels.
        for ss = 1:nSubprimaries
            Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(measureTargetSettings(ss,1,tt)*(subprimaryNInputLevels-1))); % Primary 1
            Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(measureTargetSettings(ss,2,tt)*(subprimaryNInputLevels-1))); % Primary 2
            Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(measureTargetSettings(ss,3,tt)*(subprimaryNInputLevels-1))); % Primary 3
        end
        % Measure it.
        targetSpdMeasured(:,tt) = MeasSpd(S,5,'all');
        if (options.verbose)
            fprintf('Measurement complete! - Primary %d Test Point %d ',pp,tt);
        end
    end
end

% Close PTB screen.
sca;

% Calculate LMS contrast from the measurements.
targetLMSMeasured = T_cones * targetSpdMeasured;
targetLMSContrastMeasured = ExcitationToContrast(targetLMSMeasured,bgLMS);

%% Plot the results.
% Here I used quantizedContrastImage
if (options.verbose)
    figure; hold on
    % L contrast.
    plot(measureXPixelPoints,100*measureDesiredLMSContrast(1,:),'r+','MarkerSize',4); % Measured contrast.
    plot(measureXPixelPoints,100*targetLMSContrastMeasured(1,:),'r.','MarkerSize',5); % Desired contrast.
    plot(1:imageN,contrastImage(centerN,:,1),'r','LineWidth',0.5); % Whole contrast range of the gabor patch.
    % M contrast.
    plot(measureXPixelPoints,100*measureDesiredLMSContrast(2,:),'g+','MarkerSize',4); % Measured contrast.
    plot(measureXPixelPoints,100*targetLMSContrastMeasured(2,:),'g.','MarkerSize',5); % Desired contrast.
    plot(1:imageN,contrastImage(centerN,:,2),'r','LineWidth',0.5); % Whole contrast range of the gabor patch.
    % S contrast.
    plot(measureXPixelPoints,100*measureDesiredLMSContrast(2,:),'b+','MarkerSize',4); % Measured contrast.
    plot(measureXPixelPoints,100*targetLMSContrastMeasured(2,:),'b.','MarkerSize',5); % Desired contrast.
    plot(1:imageN,contrastImage(centerN,:,3),'r','LineWidth',0.5); % Whole contrast range of the gabor patch.
    title('Image Slice, LMS Cone Contrast');
    xlabel('x position (pixels)')
    ylabel('LMS Cone Contrast (%)');
end

end
