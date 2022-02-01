function [ptCldObject,standardGaborCalObject,screenCalObj,backgroundScreenPrimaryObject] = SetupPointCloudFromGabor(...
    colorDirectionParams,rawMonochromeContrastGaborCal,screenCalObj,backgroundScreenPrimaryObject,screenPrimaryChannelObject,options)
% Set up point cloud from the gabor image.
%
% Syntax:
%    [ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,...
%    rawMonochromeContrastGaborCal,screenCalObj,backgroundChannelObject,screenPrimaryChannelObject)
%
% Description:
%    This creates the point cloud that has the all possible combinations
%    of the contrasts. This will be used to create a gabor image with a
%    desired contrast.
%
%    You can choose to use either nominal primaries or measured primaries.
%
% Inputs:
%    colorDirectionParams           - Structure with the parameters to
%                                     calculate a contrast gabor image. In
%                                     this structure, there is the number
%                                     of the target contrast gabor, which
%                                     decides the number of gabor image to
%                                     create.
%    rawMonochromeContrastGaborCal  - Monochrome contrast gabor image in a
%                                     cal format to compute the gabor image
%                                     with desired contrast and color
%                                     direction.
%    screenCalObj                   - Screen calibration object.
%    backgroundScreenPrimaryObject  - The structure with the screen
%                                     background channel primaries.
%    screenPrimaryChannelObject     - The structure with the screen priamry
%                                     channel primaries. This is only
%                                     needed if you measure the channel
%                                     primaries before making the point
%                                     cloud. If you use nominal primary,
%                                     this is not needed.
%
% Outputs:
%    ptCldObject                    - Structure with the contrasts for all
%                                     possible settings using the point
%                                     cloud method.
%    standardGaborCalObject         - Structure with the gabor contrasts
%                                     and settings in a cal format.
%    screenCalObj                   - Screen calibration object. We print
%                                     out this because the channel
%                                     primaries would be updated if we
%                                     measure the settings.
%    backgroundChannelObject        - The structure with the background
%                                     channel primaries. We want to print
%                                     out because the background
%                                     excitations would be updated if we
%                                     measure the settings.
%
% Optional key/value pairs:
%    measure                        - Default set to false. It decides if
%                                     measuring the channel primaries to
%                                     calculate the point cloud. If it sets
%                                     to true, it measures.
%    warmupTimeMinutes              - This decides the time to warm up the
%                                     projector when you decide to measure
%                                     the channel primaries. It set the
%                                     time in minutes.
%    verbose                        - Boolean. Default true. Controls
%                                     plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,gka,smo           - Wrote it.
%   01/24/22  smo                   - Made it work.
%   01/31/22  smo                   - It is possible to make multiple
%                                     target contrast gabors inside this
%                                     function.
%   02/01/22  smo                   - Added the part measuring the channel
%                                     primaries before creating point
%                                     cloud.

%% Set parameters.
arguments
    colorDirectionParams
    rawMonochromeContrastGaborCal
    screenCalObj
    backgroundScreenPrimaryObject
    screenPrimaryChannelObject
    options.measure (1,1) = false
    options.warmupTimeMinutes (1,1) = 30
    options.verbose (1,1) = true
end

%% Take some parameters out from the structure.
screenBgExcitations = backgroundScreenPrimaryObject.screenBgExcitations;
screenBgSettings    = backgroundScreenPrimaryObject.screenBgSettings;

%% Measure the channel primaries here.
%
% You can choose if you use the nominal primaries or measure it. If you use
% nominal primaries, this part will be skipped.
if (options.measure)
    % Set target screen spd which will be compared with the measured results.
    targetScreenSpd = colorDirectionParams.screenCalObj.get('P_device');
    
    % Loop and measure all primaries here.
    %
    % Open up screen and radiometer.
    [window,windowRect] = OpenPlainScreen([1 1 1]');
    OpenSpectroradiometer;
    
    % Set subprimaries to desired value and wait for them to warm up to
    % steady state.
    SetChannelSettings(screenPrimaryChannelObject.screenPrimarySettings,...
        'nInputLevels',colorDirectionParams.channelNInputLevels);
    if (options.verbose)
        fprintf('Waiting for warmup time of %d minutes ...',options.warmupTimeMinutes);
    end
    pause(60 * options.warmupTimeMinutes);
    if (options.verbose)
        fprintf('done.  Measuring.\n');
    end
    
    % Measure it.
    for pp = 1:colorDirectionParams.nPrimaries
        theScreenOnePrimarySettings = zeros(colorDirectionParams.nPrimaries,1);
        theScreenOnePrimarySettings(pp) = 1;
        targetScreenSpdMeasured(:,pp) = MeasurePlainScreenSettings(theScreenOnePrimarySettings,...
            colorDirectionParams.S, window, windowRect, 'measurementOption', options.measure, 'verbose', options.verbose);
        clear theScreenOnePrimarySettings;
    end
    
    % Plot the spd results.
    if (options.verbose)
        figure; clf;
        for pp = 1:nPrimaries
            subplot(nPrimaries,1,pp); hold on;
            plot(colorDirectionParams.wls,targetScreenSpd(:,pp),'k','LineWidth',3)
            plot(colorDirectionParams.wls,targetScreenSpdMeasured(:,pp),'r','LineWidth',2);
            xlabel('Wavelength (nm)');
            ylabel('Spectral power distribution');
            legend('Target','Measured');
            title('Comparison of raw measured and desired spds');
        end
    end
    
    % Set the primaries in the calibration to the measured results.
    %
    % It's important to also set the sensor color space, because the
    % transform between sensor/primaries is cached when we set it.
    screenCalObj.set('P_device',targetScreenSpdMeasured);
    SetSensorColorSpace(screenCalObj,colorDirectionParams.T_cones,colorDirectionParams.S);
    
    % Re-calculate the screen bg excitations.
    %
    % Not that we used to acquire the background screen settings from
    % 'ptCldScreenSettingsCheckCal' for the calculation, which was found
    % from the point cloud having the zero contrast. Here we are in the
    % pre-step making the point cloud, so the settings are loaded from the
    % screenBgSettings, which are the basically the same values.
    screenBgSpd = PrimaryToSpd(screenCalObj, SettingsToPrimary(screenCalObj,screenBgSettings));
    screenBgExcitations = colorDirectionParams.T_cones * screenBgSpd;
    
    % Update the screenBgExcitations to the struct. We will print out
    % this updated structure.
    backgroundScreenPrimaryObject.screenBgExcitations = screenBgExcitations;
end

%% Set up point cloud of contrasts for all possible settings.
%
% Make a point cloud here. It will take a while.
if (options.verbose)
    disp('Starting to make contrast point cloud...');
end
[ptCldObject.contrastPtCld, ptCldObject.ptCldSettingsCal] = ...
    SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',options.verbose);
if (options.verbose)
    disp('Contrast point cloud has been successfully created!');
end

%% Get cone contrast/excitation gabor image.
%
% Scale target cone contrast vector at max excursion by contrast modulation
% at each pixel.  This is done by a single matrix multiply plus a lead
% factor.  We work cal format here as that makes color transforms
% efficient.
nContrastPoints = size(colorDirectionParams.spatialGaborTargetContrast,2);

% Make a loop for the number of target gabor contrasts.
for cc = 1:nContrastPoints
    desiredContrastGaborCal = colorDirectionParams.spatialGaborTargetContrast(cc) * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastGaborCal;
    
    % Convert cone contrast to excitations
    desiredExcitationsGaborCal = ContrastToExcitation(desiredContrastGaborCal,screenBgExcitations);
    
    % Get primaries using standard calibration code, and desired spd without
    % quantizing.
    standardPrimariesGaborCal = SensorToPrimary(screenCalObj,desiredExcitationsGaborCal);
    desiredSpdGaborCal = PrimaryToSpd(screenCalObj,standardPrimariesGaborCal);
    
    % Gamma correct and quantize (if gamma method set to 2 above; with gamma
    % method set to zero there is no quantization).  Then convert back from
    % the gamma corrected settings.
    standardSettingsGaborCal = PrimaryToSettings(screenCalObj,standardPrimariesGaborCal);
    standardPredictedPrimariesGaborCal = SettingsToPrimary(screenCalObj,standardSettingsGaborCal);
    standardPredictedExcitationsGaborCal = PrimaryToSensor(screenCalObj,standardPredictedPrimariesGaborCal);
    standardPredictedContrastGaborCal = ExcitationsToContrast(standardPredictedExcitationsGaborCal,screenBgExcitations);
    
    % Save the results in a struct.
    standardGaborCalObject.desiredContrastGaborCal{cc} = desiredContrastGaborCal;
    standardGaborCalObject.desiredExcitationsGaborCal{cc} = desiredExcitationsGaborCal;
    standardGaborCalObject.standardPrimariesGaborCal{cc} = standardPrimariesGaborCal;
    standardGaborCalObject.desiredSpdGaborCal{cc} = desiredSpdGaborCal;
    standardGaborCalObject.standardSettingsGaborCal{cc} = standardSettingsGaborCal;
    standardGaborCalObject.standardPredictedPrimariesGaborCal{cc} = standardPredictedPrimariesGaborCal;
    standardGaborCalObject.standardPredictedExcitationsGaborCal{cc} = standardPredictedExcitationsGaborCal;
    standardGaborCalObject.standardPredictedContrastGaborCal{cc} = standardPredictedContrastGaborCal;
end
end
