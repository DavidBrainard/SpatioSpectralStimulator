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
%    lightVer                       - Deafult to true. Print out less variables
%                                     saved in the structure. It does not affect
%                                     making final gabor images, but saving some
%                                     time and memory.
%    printGaborSpds                 - Default to false. Print out the spds
%                                     of gabor images so that we can create
%                                     sRGB images easily if we want.
%    verbose                        - Boolean. Default true. Controls
%                                     plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%    01/21/22  dhb,gka,smo           - Wrote it.
%    01/24/22  smo                   - Made it work.
%    01/31/22  smo                   - It is possible to make multiple
%                                      target contrast gabors inside this
%                                      function.
%    02/01/22  smo                   - Added the part measuring the channel
%                                      primaries before creating point
%                                      cloud.
%    02/08/22  smo                   - Added an option to print out less variable
%                                      saved in the final structure.
%    05/09/22  smo                   - Added an option to make a phase shift
%                                      on sine image.
%    07/12/22  smo                   - We print out spds of gabor images in
%                                      cal format so that we can create
%                                      sRGB images easily if we want.

%% Set parameters.
arguments
    colorDirectionParams
    rawMonochromeContrastGaborCal
    screenCalObj
    backgroundScreenPrimaryObject
    screenPrimaryChannelObject
    options.measure (1,1) = false
    options.warmupTimeMinutes (1,1) = 0
    options.verbose (1,1) = true
    options.lightVer (1,1) = true
    options.printGaborSpds (1,1) = false
end

%% Take some parameters out from the structure.
screenBgExcitations = backgroundScreenPrimaryObject.screenBgExcitations;
screenBgSettings    = backgroundScreenPrimaryObject.screenBgSettings;

%% Measure the channel primaries here.
%
% You can choose if you use the nominal primaries or measure it. If you use
% nominal primaries, this part will be skipped.
if (options.measure)
    % Measurement happens here.
    %
    % As we are going to make the test images with all spatial frequencies
    % at once, we will measure screen primaries only once and save it, then
    % load it for the other cases.
    %
    % Strategy here is to save the measurement results as a separate file
    % and load it. The measurement happens if there is no such file with the
    % name.
    %
    % Get measurement file name.
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        testFilename = fullfile(testFiledir,'TestImages','MeasurementData','targetScreenSpdMeasured.mat');
    end
    
    % Load it if the file exists.
    if isfile(testFilename)
        data = load(testFilename);
        targetScreenSpdMeasured = data.targetScreenSpdMeasured;
        fprintf('Measurement file found, so skipping the measurement!');
    else
        % Measure it if the file does not exist.
        % Set target screen spd which will be compared with the measured results.
        targetScreenSpd = screenCalObj.get('P_device');
        
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
 
        % Loop and measure all primaries here.
        nPrimaries = screenCalObj.cal.nDevices;
        for pp = 1:nPrimaries
            theScreenOnePrimarySettings = zeros(nPrimaries,1);
            theScreenOnePrimarySettings(pp) = 1;
            targetScreenSpdMeasured(:,pp) = MeasurePlainScreenSettings(theScreenOnePrimarySettings,...
                colorDirectionParams.S, window, windowRect, 'measurementOption', options.measure, 'verbose', options.verbose);
            clear theScreenOnePrimarySettings;
        end
        
        % Save the measurement results.
        save(testFilename,'targetScreenSpdMeasured');
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
nPhaseShifts = length(rawMonochromeContrastGaborCal);

% Make a loop for the number of target gabor contrasts and also the number
% of phase shifts.
for ss = 1:nPhaseShifts
    rawMonochromeContrastGaborCalSingle = rawMonochromeContrastGaborCal{ss};
    
    for cc = 1:nContrastPoints
        desiredContrastGaborCal = colorDirectionParams.spatialGaborTargetContrast(cc) * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastGaborCalSingle;
        
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
        standardGaborCalObject.desiredContrastGaborCal{ss,cc} = desiredContrastGaborCal;
        standardGaborCalObject.standardSettingsGaborCal{ss,cc} = standardSettingsGaborCal;
        
        % Print out Spd if you want to make sRGB gabor image.
        if (options.printGaborSpds)
            standardGaborCalObject.desiredSpdGaborCal{ss,cc} = desiredSpdGaborCal;
        end
        
        % We will not print out these if it is run lightver.
        if (~options.lightVer)
            standardGaborCalObject.desiredExcitationsGaborCal{ss,cc} = desiredExcitationsGaborCal;
            standardGaborCalObject.standardPrimariesGaborCal{ss,cc} = standardPrimariesGaborCal;
            standardGaborCalObject.standardPredictedPrimariesGaborCal{ss,cc} = standardPredictedPrimariesGaborCal;
            standardGaborCalObject.standardPredictedExcitationsGaborCal{ss,cc} = standardPredictedExcitationsGaborCal;
            standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc} = standardPredictedContrastGaborCal;
        end
        
        fprintf('Setting up point cloud from gabor - Contrast point (%d/%d) \n', cc, nContrastPoints);
    end
    
    fprintf('Setting up point cloud from gabor - Phase shift (%d/%d) \n', ss, nPhaseShifts);
end
end
