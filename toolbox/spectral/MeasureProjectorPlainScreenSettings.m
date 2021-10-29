function [testSpdMeasured] = MeasureProjectorPlainScreenSettings(testProjectorSettings,S,options)
% Measure the SPD over projector primary settings.
%
% Syntax: [testSpdMeasured] = MeasureProjectorPlainScreenSettings(testProjectorSettings,S)
%
% Description:
%    This measures the SPD according to projector primary settings. In case
%    you use multi-subprimary device, you should set subprimary settings in
%    advance before calling this function.
%
% Inputs:
%    testProjectorSettings -      Projector input settings that reproduce
%                                 the desired contrast.
%    S -                          Spectrum measurement range in wavelength.
%
% Outputs:
%    testSpdMeasured -            Measurement results of the SPDs for the
%                                 contrast testing points.
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
%    10/19/21  smo                Now it does not take subprimary settings
%                                 as input. It only takes care of the
%                                 projector primary settings.
%    10/26/21  smo                Delete calstruct data in input, and add
%                                 'S' as input for getting spectrum
%                                 measurement range.

%% Set parameters.
arguments
    testProjectorSettings
    S
    T_cones
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

%% Followings are Disabled since 10/19/2021 (Settings including subprimary)
% Get number of discrete input levels for the device.
% So, control values go from 0 to (subprimaryNInputLevels-1).
% subprimaryNInputLevels = size(subPrimaryCalstructData{1}.get('gammaInput'),1);

% Set primary and subprimary numbers.
% nPrimaries = 3;
% nSubprimaries = subPrimaryCalstructData{1}.get('nDevices');

% Set the working range of subprimary channels (chanel 8 is not working at the moment).
% logicalToPhysical = [0:7 9:15];

% Set subprimary settings.  These ought to match what we get from the
% primary settings we computed that produced SpectralTestCal.  But
% we also ought to take the measured primaries into account, which we
% are not currently doing.
% isolatingPrimariesSpd = projectorCalObj.get('P_device');
% for pp = 1:nPrimaries
%     isolatingPrimaries(:,pp) = SpdToPrimary(subPrimaryCalstructData{pp},isolatingPrimariesSpd(:,pp));
%     subPrimarySettings(:,pp) = PrimaryToSettings(subPrimaryCalstructData{pp},isolatingPrimaries(:,pp));
% end

%% Display the projector image and measure it.
% Set projector setting for each test point on the gabor patch and measure it.
if (options.measurementOption)
    % Following subprimary setting Disabled since 10/19/2021
    % Set subprimary settings to reproduce the isolating Spd.
    %     for ss = 1:nSubprimaries
    %         Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(subPrimarySettings(ss,1)*(subprimaryNInputLevels-1))); % Primary 1
    %         Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(subPrimarySettings(ss,2)*(subprimaryNInputLevels-1))); % Primary 2
    %         Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(subPrimarySettings(ss,3)*(subprimaryNInputLevels-1))); % Primary 3
    %     end
   
    % Measure the test points.
    for tt = 1:nTestPoints
        % Set the projector settings and display it as a plane screen.
        testProjectorDisplayColor = testProjectorSettings(:,tt); % Set projector color with contrast testing point.
        OpenProjectorPlainScreen(testProjectorDisplayColor);
        % Measure it.
        testSpdMeasured(:,tt) = MeasSpd(S,5,'all');
        if (options.verbose)
            fprintf('           Measurement complete! - Test Point (%d/%d) \n',tt,nTestPoints);
        end
    end
else
    % Just print zero spectrums for both test and background
    % when skipping the measurements.
    testSpdMeasured = zeros(S(3),nTestPoints);
    if (options.verbose)
        fprintf('           Measurement has been skipped \n');
    end
end

% Close PTB screen.
CloseProjectorScreen;

%% Plot the results.
if (options.verbose)
    figure; hold on
    wls = SToWls(S); % Spectrum range.
    plot(wls,testSpdMeasured,'r-','LineWidth',2); % Test points Spd.
    title('Measured SPDs');
    xlabel('Wavelength (nm)')
    ylabel('Spectral Intensity');
end
end