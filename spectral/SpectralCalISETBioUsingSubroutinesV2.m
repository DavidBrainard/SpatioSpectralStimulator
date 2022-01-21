% SpectralCalISETBioUsingSubroutines
%
% Description:
%    Illustrate how to set up ISETBio based stimuli for the projector.
%    This version illustrates use of our encapsulated functions.
%
%    See SpectralCalCompute etc. (see "See also" list below) for a more
%    elaborated version of this the underlying computations.
%
% See also: SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%           SpectralCalISETBio

% History:
%    01/18/2022  dhb,smo    Started on it.

%% Clear
clear; close all;

%% Set key stimulus parameters.
%
% Set up color direction parameters by its condition name.
conditionName = 'LminusMSmooth';
colorDirectionParams = SetupColorDirection(conditionName);

% Set to true to get more output.
VERBOSE = false;

%% Do all calibraiton loading
[...] = LoadAndSetExperimentCalFiles(....)

%% Load screen calibration and refit its gamma.
%
% Load screen calibration.
screenCalObj = LoadCalibration(colorDirectionParams.screenCalName,...
    colorDirectionParams.screenNInputLevels,'setGammaFitMethod',true);

% Load channel calibration.
nScreenPrimaries = size(colorDirectionParams.channelCalNames,2);
for pp = 1:nScreenPrimaries
    channelCalObjs{pp} = LoadCalibration(colorDirectionParams.channelCalNames{pp},...
        colorDirectionParams.channelNInputLevels,'setGammaFitMethod',false);
end

%% Get out some data to work with.
%
% This is from the channel calibration file.
for pp = 1:nScreenPrimaries
    Scheck(pp,:) = channelCalObjs{pp}.get('S');
end
if (any(colorDirectionParams.S ~= Scheck))
    error('Mismatch between calibration file S and that specified at top');
end
wls = SToWls(colorDirectionParams.S);
nChannels = channelCalObjs{1}.get('nDevices');

%% Image spatial parameters.
%
% Image will be centered in display.
sineFreqCyclesPerDeg = 1;
gaborSdDeg = 1.5;
stimulusSizeDeg = 7;

%% Use quantized conversion from here on.
%
% Comment in the line that refits the gamma to see
% effects of extreme quantization one what follows.
%
% CalibrateFitGamma(channelCalObjs{1},10);
channelGammaMethod = 2;
for cc = 1:3
    SetGammaMethod(channelCalObjs{cc},channelGammaMethod);
end

%% Use extant machinery to get primaries from spectrum.
%
% Define wavelength range that will be used to enforce the smoothnes
% through the projection onto an underlying basis set.  We don't the whole
% visible spectrum as putting weights on the extrema where people are not
% sensitive costs us smoothness in the spectral region we care most about.
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(wls > lowProjectWl & wls < highProjectWl);

%% Find background primaries to acheive desired xy at intensity scale of display.
%
% Set parameters for getting desired background primaries.
primaryHeadRoom = 0;
targetLambda = 3;
targetBgXYZ = xyYToXYZ([colorDirectionParams.targetBgxy ; 1]);

% Adjust these to keep background in gamut
primaryBackgroundScaleFactor = 0.5;
screenBackgroundScaleFactor = 0.5;

% Make a loop for getting background for all primaries.
% Passing true for key 'Scale' causes these to be scaled reasonably
% relative to gamut, which is why we can set the target luminance
% arbitrarily to 1 just above. The scale factor determines where in the
% approximate channel gamut we aim the background at.
for pp = 1:nScreenPrimaries
    [channelBackgroundPrimaries(:,pp),channelBackgroundSpd(:,pp),channelBackgroundXYZ(:,pp)] = ...
        FindBgChannelPrimaries(targetBgXYZ, colorDirectionParams.T_xyz, channelCalObjs{pp}, colorDirectionParams.B_natural{pp}, ...
        projectIndices, primaryHeadRoom, targetLambda, 'scaleFactor', 0.6, 'Scale', true, 'Verbose', VERBOSE);
end
if (any(channelBackgroundPrimaries < 0) | any(channelBackgroundPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end
fprintf('Background primary min: %0.2f, max: %0.2f, mean: %0.2f\n', ...
    min(channelBackgroundPrimaries(:)), max(channelBackgroundPrimaries(:)), mean(channelBackgroundPrimaries(:)));

%% Find primaries with desired LMS contrast.
%
% Get isolating primaries for all screen primaries.
for pp = 1:nScreenPrimaries
    % The ambient with respect to which we compute contrast is from all
    % three primaries, which we handle via the extraAmbientSpd key-value
    % pair in the call.  The extra is for the primaries not being found in
    % the current call - the contribution from the current primary is known
    % because we pass the primaries for the background.
    otherPrimaries = setdiff(1:nScreenPrimaries,pp);
    extraAmbientSpd = 0;
    for oo = 1:length(otherPrimaries)
        extraAmbientSpd = extraAmbientSpd + channelBackgroundSpd(:,otherPrimaries(oo));
    end
    
    % Get isolating screen primaries.
    [screenPrimaryPrimaries(:,pp),screenPrimaryPrimariesQuantized(:,pp),screenPrimarySpd(:,pp),screenPrimaryContrast(:,pp),screenPrimaryModulationPrimaries(:,pp)] ...
        = FindChannelPrimaries(colorDirectionParams.targetScreenPrimaryContrastDir(:,pp), ...
        colorDirectionParams.targetPrimaryHeadroom,colorDirectionParams.targetScreenPrimaryContrasts(pp),channelBackgroundPrimaries(:,pp), ...
        colorDirectionParams.T_cones,channelCalObjs{pp},colorDirectionParams.B_natural{pp},projectIndices,colorDirectionParams.primaryHeadroom,...
        colorDirectionParams.targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
    
    % We can wonder about how close to gamut our primaries are.  Compute
    % that here.
    primaryGamutScaleFactor(pp) = MaximizeGamutContrast(screenPrimaryModulationPrimaries(:,pp),channelBackgroundPrimaries(:,pp));
    fprintf('\tPrimary %d, gamut scale factor is %0.3f\n',pp,primaryGamutScaleFactor(pp));
    
    % Find the channel settings that correspond to the desired screen
    % primaries.
    screenPrimarySettings(:,pp) = PrimaryToSettings(channelCalObjs{pp},screenPrimaryPrimaries(:,pp));
end

%% How close are spectra to subspace defined by basis?
%
% This part has been updated using the loop to make it short.
for pp = 1:nScreenPrimaries
    isolatingNaturalApproxSpd(:,pp) = colorDirectionParams.B_natural{pp} * (colorDirectionParams.B_natural{pp}(projectIndices,:)\screenPrimarySpd(projectIndices,pp));
end

% Plot of the screen primary spectra.
figure; clf;
for pp = 1:nScreenPrimaries
    subplot(2,2,pp); hold on;
    plot(wls,screenPrimarySpd(:,pp),'b','LineWidth',2);
    plot(wls,isolatingNaturalApproxSpd(:,pp),'r:','LineWidth',1);
    plot(wls(projectIndices),screenPrimarySpd(projectIndices,pp),'b','LineWidth',4);
    plot(wls(projectIndices),isolatingNaturalApproxSpd(projectIndices,pp),'r:','LineWidth',3);
    xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
    title(append('Primary ', num2str(pp)));
end

%% Set the screen primaries.
%
% We want these to match those we set up with the channel calculations
% above.  Need to reset sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device',screenPrimarySpd);
SetSensorColorSpace(screenCalObj,colorDirectionParams.T_cones,colorDirectionParams.S);

%% Set screen gamma method.
%
% If we set to 0, there is no quantization and the result is excellent.
% If we set to 2, this is quantized at 256 levels and the result is more
% of a mess.  The choice of 2 represents what we think will actually happen
% since the real device is quantized.
%
% The point cloud method below reduces this problem.
screenGammaMethod = 2;
SetGammaMethod(screenCalObj,screenGammaMethod);

%% Create ISETBio display from the calibration file.
%
% On the DMD size, Derek writes (email 2022-01-11):
%     Looking at the system design, I find that the conversion is 7.74350
%     degrees : 12 mm = 0.64529 degrees per mm. This indicates that the
%     diagonal of the DMD occupies a FOV of 2*7.74350 = 15.487 degrees ~
%     15.5 degrees.
%
%     We can use this to solve for the h and v size in degrees.
%     Let the horizontal size by x and the diagonal be d.  We know
%     that d^2 = x^2 * (1+(screenVertPixels/screenHorizPixels)^2). So
%     x = sqrt(d^2/(1+(screenVertPixels/screenHorizPixels)^2)).
%
% The DMD dimensions and distance are dummied up so that the visual
% angle matches that of our optical system, but so that the optical
% distance is large enough to mimic optical infinity.
%
% Display parameters.
screenDiagSizeDeg = 15.5;
screenDistanceVirtualMeters = 10;
screenSizePixels = [1920 1080];
screenDiagSizePixels = vecnorm(screenSizePixels);

% Convert the screen size from degrees to meters/inches.
%
% Note that we need to do this along the diagonal because degrees aren't
% linear in meters, so we want to work first in the physical units of the
% display, not in degrees.
[screenDiagSizeMeters,screenSizeMeters,screenSizeInches] = ...
    DegToMeters(screenDiagSizeDeg,screenDistanceVirtualMeters,screenSizePixels,'DegToInches',true);

% Get horizontal and vertical size of screen in degrees. We take pixels per
% degree along the diagonal as the best compromise, need to use that when
% we compute image sizes below.
screenSizeDeg = 2 * atand(screenSizeMeters/(2*screenDistanceVirtualMeters));
screenPixelsPerDeg = screenDiagSizePixels / screenDiagSizeDeg;

% Get dpi and make sure everything is consistent.
screenDpi = vecnorm(screenSizePixels)/vecnorm(screenSizeInches);
screenDpiChk = mean(screenSizePixels ./ screenSizeInches);
if (abs((screenDpi - vecnorm(screenDpiChk))/screenDpi) > 1e-6)
    error('Screen is not rigid');
end
inchesPerMeter = 39.3701;
screenDpm = screenDpi*inchesPerMeter;

% Create ISETBio display object here.
%
% Note that the function takes T_cones, S, screenGammaMethod for checking
% if the reverse calculation (from ISETBio to calibration object) can get
% the same values.
[ISETBioDisplayObject,screenCalObjFromISETBio] = MakeISETBioDisplayObj(screenCalObj,screenDistanceVirtualMeters,...
    screenSizeMeters,screenSizePixels,colorDirectionParams.T_cones,colorDirectionParams.S,screenGammaMethod,'verbose',VERBOSE);

%% Set up desired background.
%
% We aim for the background that we said we wanted when we built the screen primaries.
desiredBgExcitations = screenBackgroundScaleFactor * colorDirectionParams.T_cones * sum(channelBackgroundSpd,2);
screenBgSettings = SensorToSettings(screenCalObj,desiredBgExcitations);
screenBgExcitations = SettingsToSensor(screenCalObj,screenBgSettings);

% Plot it.
figure; clf; hold on;
plot(desiredBgExcitations,screenBgExcitations,'ro','MarkerFaceColor','r','MarkerSize',12);
axis('square');
xlim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
ylim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
xlabel('Desired bg excitations'); ylabel('Obtained bg excitations');
title('Check that we obtrain desired background excitations');
fprintf('Screen settings to obtain background: %0.2f, %0.2f, %0.2f\n', ...
    screenBgSettings(1),screenBgSettings(2),screenBgSettings(3));

%% Make a monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image. Multiply by the max
% contrast vector to get the LMS contrast image.
fprintf('Making Gabor contrast image\n');
[rawMonochromeUnquantizedContrastGaborImage, rawMonochromeUnquantizedContrastGaborCal, ...
    stimulusN, centerN, stimulusHorizSizeDeg, stimulusHorizSizeMeters] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenPixelsPerDeg,screenDpm,'verbose',VERBOSE);

%% Quantize the contrast image to a (large) fixed number of levels.
%
% This allows us to speed up the image conversion without any meaningful
% loss of precision. If you don't like it, increase number of quantization
% bits until you are happy again.
nQuantizeBits = 14;
nQuantizeLevels = 2^nQuantizeBits;
rawMonochromeContrastGaborCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastGaborCal+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;

% Plot of how quantized version does in obtaining desired contrats.
figure; clf;
plot(rawMonochromeUnquantizedContrastGaborCal(:),rawMonochromeContrastGaborCal(:),'r+');
axis('square');
xlim([0 1]); ylim([0 1]);
xlabel('Unquantized Gabor contrasts');
ylabel('Quantized Gabor contrasts');
title('Effect of contrast quantization');

%% Get cone contrast/excitation gabor image.
%
% Scale target cone contrast vector at max excursion by contrast modulation
% at each pixel.  This is done by a single matrix multiply plus a lead
% factor.  We work cal format here as that makes color transforms
% efficient.
desiredContrastGaborCal = colorDirectionParams.spatialGaborTargetContrast * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastGaborCal;

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

%% Set up point cloud of contrasts for all possible settings
[contrastPtCld, ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',VERBOSE);

%% Get image from point cloud in cal format.
uniqueQuantizedSettingsGaborCal = SettingsFromPointCloud(contrastPtCld,desiredContrastGaborCal,ptCldSettingsCal);

% Print out min/max of settings
fprintf('Gabor image min/max settings: %0.3f, %0.3f\n',min(uniqueQuantizedSettingsGaborCal(:)), max(uniqueQuantizedSettingsGaborCal(:)));

% Get contrasts we think we have obtianed
uniqueQuantizedExcitationsGaborCal = SettingsToSensor(screenCalObj,uniqueQuantizedSettingsGaborCal);
uniqueQuantizedContrastGaborCal = ExcitationsToContrast(uniqueQuantizedExcitationsGaborCal,screenBgExcitations);

% Plot of how well point cloud method does in obtaining desired contrats
figure; clf;
plot(desiredContrastGaborCal(:),uniqueQuantizedContrastGaborCal(:),'r+');
axis('square');
xlabel('Desired L, M or S contrast');
ylabel('Predicted L, M, or S contrast');
title('Quantized unique point cloud image method');

%% Convert representations we want to take forward to image format
desiredContrastGaborImage = CalFormatToImage(desiredContrastGaborCal,stimulusN,stimulusN);
standardPredictedContrastImage = CalFormatToImage(standardPredictedContrastGaborCal,stimulusN,stimulusN);
standardSettingsGaborImage = CalFormatToImage(standardSettingsGaborCal,stimulusN,stimulusN);
uniqueQuantizedContrastGaborImage = CalFormatToImage(uniqueQuantizedContrastGaborCal,stimulusN,stimulusN);

%% Put the image into an ISETBio scene
%
% These calls are a bit slow for large images and the fine wavelength
% sampling used here. But these would be done as pre-compute steps so
% it doesn't seem worth trying to optimize at this point.
ISETBioGaborScene = sceneFromFile(standardSettingsGaborImage,'rgb', [], ISETBioDisplayObject);
sceneWindow(ISETBioGaborScene);

% Check stimulus dimensions match. These are good to about a percent, which
% we can live with.
stimulusHorizSizeMetersChk = sceneGet(ISETBioGaborScene,'width');
stimulusHorizSizeDegChk = sceneGet(ISETBioGaborScene,'horizontal fov');
if (abs(stimulusHorizSizeMeters - stimulusHorizSizeMetersChk)/stimulusHorizSizeMeters > 0.01)
    error('Horizontal size in meters mismatch of too much');
end
if (abs(stimulusHorizSizeDeg - stimulusHorizSizeDegChk)/stimulusHorizSizeDeg > 0.01)
    error('Horizontal size in deg mismatch of too much');
end

% Calculate cone excitations from the ISETBio scene.
% These should match what we get when we compute
% outside of ISETBio. And indeed!

% ISETBio energy comes back as power per nm, we need to convert to power
% per wlband to work with PTB, by multiplying by S(2).
ISETBioGaborImage = sceneGet(ISETBioGaborScene,'energy') * colorDirectionParams.S(2);
[ISETBioGaborCal,ISETBioM,ISETBioN] = ImageToCalFormat(ISETBioGaborImage);
ISETBioPredictedExcitationsGaborCal = colorDirectionParams.T_cones * ISETBioGaborCal;
limMin = 0.01; limMax = 0.02;
figure; clf; hold on;
plot(standardPredictedExcitationsGaborCal(1,:), ISETBioPredictedExcitationsGaborCal(1,:),'r+');
plot(standardPredictedExcitationsGaborCal(2,:), ISETBioPredictedExcitationsGaborCal(2,:),'g+');
plot(standardPredictedExcitationsGaborCal(3,:), ISETBioPredictedExcitationsGaborCal(3,:),'b+');
plot([limMin limMax], [limMin limMax]);
xlabel('Standard Cone Excitations');
ylabel('ISETBio Cone Excitations');
axis('square'); xlim([limMin limMax]); ylim([limMin limMax]);
title('Cone Excitations Comparison');
if (max(abs(standardPredictedExcitationsGaborCal(:) - ISETBioPredictedExcitationsGaborCal(:)) ./ ...
        standardPredictedExcitationsGaborCal(:)) > 1e-6)
    error('Standard and ISETBio data do not agree well enough');
end

% Go back to the RGB image starting with the ISETBio representation.
primaryFromISETBioGaborCal = screenCalObjFromISETBio.get('P_device') \ ...
    (ISETBioGaborCal-screenCalObjFromISETBio.get('P_ambient'));
settingsFromISETBioGaborCal = PrimaryToSettings(screenCalObjFromISETBio,primaryFromISETBioGaborCal);
if (max(abs(standardSettingsGaborCal(:)-settingsFromISETBioGaborCal(:))./standardSettingsGaborCal(:)) > 1e-6)
    error('Cannot get home again in settings land');
end

%% SRGB image via XYZ, scaled to display
predictedXYZCal = colorDirectionParams.T_xyz * desiredSpdGaborCal;
SRGBPrimaryCal = XYZToSRGBPrimary(predictedXYZCal);
scaleFactor = max(SRGBPrimaryCal(:));
SRGBCal = SRGBGammaCorrect(SRGBPrimaryCal/(2*scaleFactor),0);
SRGBImage = uint8(CalFormatToImage(SRGBCal,stimulusN,stimulusN));

% Show the SRGB image
figure; imshow(SRGBImage);
title('SRGB Gabor Image');

%% Show the settings image
figure; clf;
imshow(standardSettingsGaborImage);
title('Image of settings');

%% Plot slice through predicted LMS contrast image.
%
% Note that the y-axis in this plot is individual cone contrast, which is
% not the same as the vector length contrast of the modulation.
plotAxisLimit = 100 * colorDirectionParams.spatialGaborTargetContrast;

figure; hold on
plot(1:stimulusN,100*standardPredictedContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:stimulusN,100*standardPredictedContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:stimulusN,100*standardPredictedContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
if (screenGammaMethod == 2)
    title('Image Slice, SensorToSettings Method, Quantized Gamma, LMS Cone Contrast');
else
    title('Image Slice, SensorToSettings Method, No Quantization, LMS Cone Contrast');
end
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Plot slice through point cloud LMS contrast image.
%
% Note that the y-axis in this plot is individual cone contrast, which is
% not the same as the vector length contrast of the modulation.
figure; hold on
plot(1:stimulusN,100*uniqueQuantizedContrastGaborImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:stimulusN,100*uniqueQuantizedContrastGaborImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:stimulusN,100*uniqueQuantizedContrastGaborImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:stimulusN,100*desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, Point Cloud Method, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Generate some settings values corresponding to known contrasts
%
% The reason for this is to measure and check these.  This logic follows
% how we handled an actual gabor image above. We don't actually need to
% quantize to 14 bits here on the contrast, but nor does it hurt.
rawMonochromeUnquantizedContrastCheckCal = [0 0.25 -0.25 0.5 -0.5 1 -1];
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
desiredContrastCheckCal = colorDirectionParams.spatialGaborTargetContrast * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastCheckCal;
desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,screenBgExcitations);

% For each check calibration find the settings that
% come as close as possible to producing the desired excitations.
%
% If we measure for a uniform field the spectra corresopnding to each of
% the settings in the columns of ptCldScreenSettingsCheckCall, then
% compute the cone contrasts with respect to the backgound (0 contrast
% measurement, first settings), we should approximate the cone contrasts in
% desiredContrastCheckCal.
ptCldScreenSettingsCheckCal = SettingsFromPointCloud(contrastPtCld,desiredContrastCheckCal,ptCldSettingsCal);
ptCldScreenPrimariesCheckCal = SettingsToPrimary(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenSpdCheckCal = PrimaryToSpd(screenCalObj,ptCldScreenPrimariesCheckCal);
ptCldScreenExcitationsCheckCal = SettingsToSensor(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenContrastCheckCal = ExcitationsToContrast(ptCldScreenExcitationsCheckCal,screenBgExcitations);
figure; clf; hold on;
plot(desiredContrastCheckCal(:),ptCldScreenContrastCheckCal(:),'ro','MarkerSize',10,'MarkerFaceColor','r');
xlim([0 plotAxisLimit/100]); ylim([0 plotAxisLimit/100]); axis('square');
xlabel('Desired'); ylabel('Obtained');
title('Check of desired versus obtained check contrasts');

% Check that we can recover the settings from the spectral power
% distributions, etc.  This won't necessarily work perfectly, but should be
% OK.
for tt = 1:size(ptCldScreenSettingsCheckCal,2)
    ptCldPrimariesFromSpdCheckCal(:,tt) = SpdToPrimary(screenCalObj,ptCldScreenSpdCheckCal(:,tt),'lambda',0);
    ptCldSettingsFromSpdCheckCal(:,tt) = PrimaryToSettings(screenCalObj,ptCldScreenSettingsCheckCal(:,tt));
end
figure; clf; hold on
plot(ptCldScreenSettingsCheckCal(:),ptCldSettingsFromSpdCheckCal(:),'+','MarkerSize',12);
xlim([0 1]); ylim([0 1]);
xlabel('Computed primaries'); ylabel('Check primaries from spd'); axis('square');

% Make sure that screenPrimarySettings leads to screenPrimarySpd
clear screenPrimarySpdCheck
for pp = 1:length(channelCalObjs)
    screenPrimarySpdCheck(:,pp) = PrimaryToSpd(channelCalObjs{pp},SettingsToPrimary(channelCalObjs{pp},screenPrimarySettings(:,pp)));
end
figure; clf; hold on
plot(wls,screenPrimarySpdCheck,'k','LineWidth',4);
plot(wls,screenPrimarySpd,'r','LineWidth',2);
xlabel('Wavelength'); ylabel('Radiance');
title('Check of consistency between screen primaries and screen primary spds');

%% Save out what we need to check things on the DLP
screenSettingsImage = standardSettingsGaborImage;
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    save(testFilename,'S','T_cones','screenCalObj','channelCalObjs','screenSettingsImage', ...
        'screenPrimaryPrimaries','screenPrimarySettings','screenPrimarySpd',...
        'desiredContrastCheckCal', ...
        'ptCldScreenSettingsCheckCal','ptCldScreenContrastCheckCal','ptCldScreenSpdCheckCal', ...
        'nQuantizeLevels','screenNInputLevels','targetStimulusContrastDir','spatialGaborTargetContrast');
end
