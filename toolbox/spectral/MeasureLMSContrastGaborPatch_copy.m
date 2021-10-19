function [testSpdMeasured,bgSpdMeasured] = MeasureLMSContrastGaborPatch_copy(testProjectorSettings,bgProjectorSettings,...
    projectorCalObj,subPrimaryCalstructData,T_cones,subprimaryNInputLevels,options)
% Measure the LMS contrasts at some points on the gabor patch image.
%
% Syntax: [testSpdMeasured] = MeasureLMSContrastGaborPatch_copy(testProjectorSettings,bgProjectorSettings,...
%                                                               projectorCalObj,subPrimaryCalstructData,T_cones,subprimaryNInputLevels)
%
% Description:
%    This measures and calculates the SPD for calculating LMS contrasts
%    of the modulated gabor patch, so we can check if the modulation
%    was done properly.
%
% Inputs:
%    testProjectorSettings -      Projector input settings that reproduce
%                                 the desired contrast.
%    bgProjectorSettings -        Projector input settings for the
%                                 background.
%    projectorCalObj -            The calibration object that contains the
%                                 the data for the projector.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with.
%    T_cones -                    Spectral cone sensitivities in standard PTB format.
%    subprimaryNInputLevels -     Device max input levels for the
%                                 subprimaries.
%
% Outputs:
%    testSpdMeasured -            Measurement results of the SPDs for the
%                                 contrast testing points.
%    bgSpdMeasured -              Measurement result of the SPD for the
%                                 bakground.
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
%    10/15/21  smo                Fixed the error, clean it, better naming
%                                 the variables.
%    10/18/21  smo                Made it simpler and it takes projector
%                                 settings and save out the spd data only.

%% Set parameters.
arguments
    testProjectorSettings
    bgProjectorSettings
    projectorCalObj
    subPrimaryCalstructData
    T_cones
    subprimaryNInputLevels
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

% Number of contrast test points.
nTestPoints = size(testProjectorSettings,2);

% Set the working range of subprimary channels (chanel 8 is not working at the moment).
logicalToPhysical = [0:7 9:15];

% Measurement range.
S = subPrimaryCalstructData{1}.get('S');

% Set primary and subprimary numbers.
nPrimaries = 3;
nSubprimaries = subPrimaryCalstructData{1}.get('nDevices');

% Set subprimary settings.
isolatingPrimariesSpd = projectorCalObj.get('P_device');
for pp = 1:nPrimaries
    isolatingPrimaries(:,pp) = SpdToPrimary(subPrimaryCalstructData{pp},isolatingPrimariesSpd(:,pp));
    subPrimarySettings(:,pp) = PrimaryToSettings(subPrimaryCalstructData{pp},isolatingPrimaries(:,pp));
end

%% Display the projector image and measure it.
% Set projector setting for each test point on the gabor patch and measure it.
if (options.measurementOption)
    % Set subprimary settings to reproduce the isolating Spd.
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(subPrimarySettings(ss,1)*(subprimaryNInputLevels-1))); % Primary 1
        Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(subPrimarySettings(ss,2)*(subprimaryNInputLevels-1))); % Primary 2
        Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(subPrimarySettings(ss,3)*(subprimaryNInputLevels-1))); % Primary 3
    end
    
    % Get ready to use PTB.
    PsychDefaultSetup(2); % PTB pre-setup
    screens = Screen('Screens');
    screenNumber = max(screens);
    
    % Measure the test points.
    for tt = 1:nTestPoints
        % Set the projector settings and display it as a plane screen.
        projectorDisplayColor = testProjectorSettings(:,tt); % Set projector color with contrast testing point.
        [window, windowRect] = PsychImaging('OpenWindow', screenNumber, projectorDisplayColor);
        % Measure it.
        testSpdMeasured(:,tt) = MeasSpd(S,5,'all');
        if (options.verbose)
            fprintf('           Measurement complete! - Test Point (%d/%d) \n',tt,nTestPoints);
        end
    end
    
    % Measure the background.
    % Set the projector settings and display it as a plane screen.
    projectorDisplayColor = bgProjectorSettings; % This part sets RGB values of the projector image.
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, projectorDisplayColor);
    % Measure it.
    bgSpdMeasured = MeasSpd(S,5,'all');
    if (options.verbose)
        fprintf('           Measurement complete! - Background \n');
    end
    
else
    % Just print zero spectrums for both test and background
    % when skipping the measurements.
    testSpdMeasured = zeros(S(3),nTestPoints);
    bgSpdMeasured = zeros(S(3),1);
    if (options.verbose)
        fprintf('           Measurement has been skipped \n');
    end
end

% Close PTB screen.
sca;

%% Plot the results.
if (options.verbose)
    figure; hold on
    wls = SToWls(S); % Spectrum range.
    plot(wls,testSpdMeasured); % Test points Spd.
    plot(wls,bgSpdMeasured,'k-'); % Background Spd.
    title('Measured SPDs');
    xlabel('Wavelength (nm)')
    ylabel('Spectral Intensity');
end
end