function [testLMSContrastMeasured,testSpdMeasured] = MeasureLMSContrastGaborPatch(contrastImage,isolatingPrimaries,projectorCalObj,bgSpd,T_cones,...
    subprimaryNInputLevels,subPrimaryCalstructData,options)
% Measure the LMS contrasts at some points on the gabor patch image.
%
% Syntax: [targetLMSContrastMeasured] = MeasureLMSContrastGaborPatch(quantizedContrastImage,isolatingPrimaries,projectorCalObj,bgSpd,...
%    subprimaryNInputLevels,subPrimaryCalstructData)
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
%    bgSpd -                      Background Spd of the gabor patch.
%    T_cones -                    Spectral cone sensitivities in standard PTB format.
%    subprimaryNInputLevels -     Device max input levels for the
%                                 subprimaries.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with.
%
% Outputs:
%    testLMSContrastMeasured -    Measurement results of the some points on
%                                 the contrast gabor patch image.
%    testSpdMeasured -            Measurement results of the SPDs for the
%                                 points on the slice of the gabor patch.
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
    T_cones
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
S = subPrimaryCalstructData{1}.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalstructData{1}.get('nDevices');

% Get background LMS.
bgLMS = T_cones * sum(bgSpd,2);

% Get test image size and set the center point.
imageN = size(contrastImage,1);
centerN = imageN/2;

%% Pick the points on the Gabor patch and calculate target settings for measurements.
% Set measurement points on the Gabor patch. Points are decided on the
% X-pixel 2-D plane.
% Now it picks five points of 30%, 40%, 50%, 60%, 70% positions of X pixel
% of the contrast image which should be changed later on.
measureXPixelPoints = round([0.3:0.1:0.7] * imageN);
nTestPoints = size(measureXPixelPoints,2);

% Calculate LMS contrast.
SetGammaMethod(projectorCalObj,0);
testLMSContrast = squeeze(contrastImage(centerN,measureXPixelPoints,:))';
testLMS = ContrastToExcitation(testLMSContrast,bgLMS);
testPrimariesRatio = SensorToPrimary(projectorCalObj,testLMS);

% Get target settings for the projector.
for tt = 1:nTestPoints
    % Set target primaries.
    testSubprimaries(:,:,tt) = isolatingPrimaries .* testPrimariesRatio(:,tt)';
    for pp = 1:nPrimaries
        % Set target settings per each test point.
        testSubprimarySettings(:,pp,tt) = PrimaryToSettings(subPrimaryCalstructData{pp},testSubprimaries(:,pp,tt));
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
            Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(testSubprimarySettings(ss,1,tt)*(subprimaryNInputLevels-1))); % Primary 1
            Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(testSubprimarySettings(ss,2,tt)*(subprimaryNInputLevels-1))); % Primary 2
            Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(testSubprimarySettings(ss,3,tt)*(subprimaryNInputLevels-1))); % Primary 3
        end
        % Measure it.
        testSpdMeasured(:,tt) = MeasSpd(S,5,'all');
        if (options.verbose)
            fprintf('           Measurement complete! - Test Point (%d/%d) \n',tt,nTestPoints);
        end
    end
else
     % Show the same results as background when measurement skipped.
       testSpdMeasured = ones(S(3),nTestPoints).* sum(bgSpd,2); 
end

% Close PTB screen.
sca;

% Calculate LMS contrast from the measurements.
testLMSMeasured = T_cones * testSpdMeasured;
testLMSContrastMeasured = ExcitationToContrast(testLMSMeasured,bgLMS);

%% Plot the results.
if (options.verbose)
    figure; hold on
    % L contrast.
    plot(measureXPixelPoints,100*testLMSContrast(1,:),'r+-','MarkerSize',4,'LineWidth',0.5); % Measured contrast.
    plot(measureXPixelPoints,100*testLMSContrastMeasured(1,:),'ro--','MarkerSize',7,'LineWidth',1); % Desired contrast.
    % M contrast.
    plot(measureXPixelPoints,100*testLMSContrast(2,:),'g+-','MarkerSize',4,'LineWidth',0.5); % Measured contrast.
    plot(measureXPixelPoints,100*testLMSContrastMeasured(2,:),'go--','MarkerSize',7,'LineWidth',1); % Desired contrast.
    % S contrast.
    plot(measureXPixelPoints,100*testLMSContrast(3,:),'b+-','MarkerSize',4,'LineWidth',0.5); % Measured contrast.
    plot(measureXPixelPoints,100*testLMSContrastMeasured(3,:),'bo--','MarkerSize',7,'LineWidth',1); % Desired contrast.
    title('Image Slice, LMS Cone Contrast');
    xlabel('x position (pixels)')
    ylabel('LMS Cone Contrast (%)');
    legend('Desired-L','Measured-L','Desired-M','Measured-M','Desired-S','Measured-S');
end
end