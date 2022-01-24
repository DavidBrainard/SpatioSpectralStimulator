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

%% Do all calibraiton loading.
screenGammaMethod = 2;
[screenCalObj,channelCalObjs] = LoadAndSetExperimentCalFiles(colorDirectionParams,'screenGammaMethod',screenGammaMethod,'verbose',VERBOSE);

%% Image spatial parameters.
%
% Image will be centered in display.
sineFreqCyclesPerDeg = 1;
gaborSdDeg = 1.5;
stimulusSizeDeg = 7;

%% Use extant machinery to get primaries from spectrum.
%
% Define wavelength range that will be used to enforce the smoothnes
% through the projection onto an underlying basis set.  We don't the whole
% visible spectrum as putting weights on the extrema where people are not
% sensitive costs us smoothness in the spectral region we care most about.
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(colorDirectionParams.wls > lowProjectWl & colorDirectionParams.wls < highProjectWl);

%% Find primaries with desired LMS contrast.
[screenPrimaryChannelObj,bgChannelObject] = SetupChannelPrimaries(colorDirectionParams,channelCalObjs,projectIndices);

%% Set the screen primaries.
%
% We want these to match those we set up with the channel calculations
% above.  Need to reset sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device',screenPrimaryChannelObj.screenPrimarySpd);
SetSensorColorSpace(screenCalObj,colorDirectionParams.T_cones,colorDirectionParams.S);

%% Create ISETBio display from the calibration file.
[ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj);

%% Background
[bgScreenPrimaryObject] = SetupBackground(colorDirectionParams,screenCalObj,bgChannelObject);

%% Make a monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image. Multiply by the max
% contrast vector to get the LMS contrast image.
[rawMonochromeUnquantizedContrastGaborImage, rawMonochromeUnquantizedContrastGaborCal, ...
    stimulusN, centerN, stimulusHorizSizeDeg, stimulusHorizSizeMeters] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,'verbose',VERBOSE);

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
[ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,screenCalObj,bgScreenPrimaryObject.screenBgExcitations);

%% Make image from point cloud
gaborImageObject = MakeImageSettingsFromPtCld(ptCldObject,screenCalObj,standardGaborCalObject,bgScreenPrimaryObject.screenBgExcitations,stimulusN);

%% Put the image into an ISETBio scene
ISETBioGaborCalObject = MakeISETBioSceneFromImage(colorDirectionParams,gaborImageObject,standardGaborCalObject,ISETBioDisplayObject,screenSizeObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg);

% Go back to the RGB image starting with the ISETBio representation.
[primaryFromISETBioGaborCal, settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject);

%% SRGB image via XYZ, scaled to display
predictedXYZCal = colorDirectionParams.T_xyz * standardGaborCalObject.desiredSpdGaborCal;
SRGBPrimaryCal = XYZToSRGBPrimary(predictedXYZCal);
scaleFactor = max(SRGBPrimaryCal(:));
SRGBCal = SRGBGammaCorrect(SRGBPrimaryCal/(2*scaleFactor),0);
SRGBImage = uint8(CalFormatToImage(SRGBCal,stimulusN,stimulusN));

% Show the SRGB image
figure; imshow(SRGBImage);
title('SRGB Gabor Image');

%% Show the settings image
figure; clf;
imshow(gaborImageObject.standardSettingsGaborImage);
title('Image of settings');

%% Plot slice through predicted LMS contrast image.
%
% Note that the y-axis in this plot is individual cone contrast, which is
% not the same as the vector length contrast of the modulation.
plotAxisLimit = 100 * colorDirectionParams.spatialGaborTargetContrast;

figure; hold on
plot(1:stimulusN,100 * gaborImageObject.standardPredictedContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:stimulusN,100 * gaborImageObject.standardPredictedContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:stimulusN,100 * gaborImageObject.standardPredictedContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
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
plot(1:stimulusN,100 * gaborImageObject.uniqueQuantizedContrastGaborImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:stimulusN,100 * gaborImageObject.uniqueQuantizedContrastGaborImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:stimulusN,100 * gaborImageObject.uniqueQuantizedContrastGaborImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:stimulusN,100 * gaborImageObject.desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
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
desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,bgScreenPrimaryObject.screenBgExcitations);

% For each check calibration find the settings that
% come as close as possible to producing the desired excitations.
%
% If we measure for a uniform field the spectra corresopnding to each of
% the settings in the columns of ptCldScreenSettingsCheckCall, then
% compute the cone contrasts with respect to the backgound (0 contrast
% measurement, first settings), we should approximate the cone contrasts in
% desiredContrastCheckCal.
ptCldScreenSettingsCheckCal = SettingsFromPointCloud(ptCldObject.contrastPtCld,desiredContrastCheckCal,ptCldObject.ptCldSettingsCal);
ptCldScreenPrimariesCheckCal = SettingsToPrimary(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenSpdCheckCal = PrimaryToSpd(screenCalObj,ptCldScreenPrimariesCheckCal);
ptCldScreenExcitationsCheckCal = SettingsToSensor(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenContrastCheckCal = ExcitationsToContrast(ptCldScreenExcitationsCheckCal, bgScreenPrimaryObject.screenBgExcitations);
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
    screenPrimarySpdCheck(:,pp) = PrimaryToSpd(channelCalObjs{pp},SettingsToPrimary(channelCalObjs{pp}, screenPrimaryChannelObj.screenPrimarySettings(:,pp)));
end
figure; clf; hold on
plot(colorDirectionParams.wls, screenPrimarySpdCheck,'k','LineWidth',4);
plot(colorDirectionParams.wls, screenPrimaryChannelObj.screenPrimarySpd,'r','LineWidth',2);
xlabel('Wavelength'); ylabel('Radiance');
title('Check of consistency between screen primaries and screen primary spds');

%% Save out what we need to check things on the DLP
screenSettingsImage = standardSettingsGaborImage;
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    save(testFilename,'colorDirectionParams');
end
