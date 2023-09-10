%{

%% Run a full set of examples
clear;
conditionNameList = {'MelDirected1'}; % 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2 1 2 5 10];
gaborSdDeg = 100;
stimulusSizeDeg = 4;

for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        lmsmelContrastNominal(:,ss,cc) = SpectralCalISETBioMel(conditionNameList{cc},sineFreqCyclesPerDegList(ss), ...
            gaborSdDeg,stimulusSizeDeg);
    end
end

%}

%{

%% Baseline example
clear;
conditionName = 'MelDirected1';
sineFreqCyclesPerDeg = 2;
gaborSdDeg = 100;
stimulusSizeDeg = 4;

lmsmelContrastNominal = SpectralCalISETBioMel(conditionName,sineFreqCyclesPerDeg, ...
    gaborSdDeg,stimulusSizeDeg);
   
%}

function lmsmelContrastNominal = SpectralCalISETBioMel(conditionName,sineFreqCyclesPerDeg,gaborSdDeg,stimulusSizeDeg)
% SpectralCalISETBioMel
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
%    01/18/22  dhb,smo    Started on it.
%    01/26/22  smo        It is working well with the substituting
%                         functions!
%    04/15/22  smo        Corrections because of the sub routines now
%                         save the variables in cell format.
%    10/18/22  smo        Updated it to use to check calibration every
%                         time we make test images for the experiment.
%    09/06/23  dhb        Melanopsin version

%% Clear.
close all;

%% Set key stimulus parameters.
%
% Set up color direction parameters by its condition name.
switch (conditionName)    
    case 'MelDirected1'
        targetScreenPrimaryContrasts = 0.15;
        spatialGaborTargetContrast = 0.20;
    case 'IsochromaticControl'
        targetScreenPrimaryContrasts = 0.75;
        spatialGaborTargetContrast = 0.75;
end

% Gamma method
%   0 -> no quantization modeled
%   2 -> with quantization
screenGammaMethod = 2;

% Call routine to get basic parameters structure set up.
coneFieldSizeDeg = 10;  % Eccentricity times 2
colorDirectionParams = SetupColorDirectionMel(conditionName,...
    'targetScreenPrimaryContrasts',targetScreenPrimaryContrasts,'spatialGaborTargetContrast',spatialGaborTargetContrast, ...
    'fieldSizeDegs', coneFieldSizeDeg);

% Set to true to get more output.
VERBOSE = true;

%% Do all calibration loading.
[screenCalObj,channelCalObjs] = LoadAndSetExperimentCalFiles(colorDirectionParams,'screenGammaMethod',screenGammaMethod,'verbose',VERBOSE);

%% Calculation quantization.  
%
% This is a calculation quantization, not a model of display
% quantization.  Making it smaller sppeds thigns up, but reduces
% accuracy of computed image.
%
% You can drop to 9 without much change and speed things up.
nQuantizeBits = 11;
nQuantizeLevels = 2^nQuantizeBits;

%% Output names and places
projectFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
sceneOutputStr = sprintf('%s_Size_%0.1f_Sf_%0.1f_Sd_%0.1f_GammaMethod_%d',conditionName,stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenGammaMethod);
outputDir = fullfile(projectFiledir,sceneOutputStr);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
sceneOutputFilename = fullfile(outputDir,'sceneOutput.mat');

%% Use extant machinery to get primaries from spectrum.
%
% Define wavelength range that will be used to enforce the smoothness
% through the projection onto an underlying basis set.  We don't the whole
% visible spectrum as putting weights on the extrema where people are not
% sensitive costs us smoothness in the spectral region we care most about.
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(colorDirectionParams.wls > lowProjectWl & colorDirectionParams.wls < highProjectWl);

%% Find primaries with desired LMS contrast.
[screenPrimaryChannelObject,backgroundChannelObject] = SetupChannelPrimariesMel(colorDirectionParams,channelCalObjs,projectIndices,'verbose',VERBOSE);

%% Set the screen primaries.
%
% We want these to match those we set up with the channel calculations
% above.  Need to reset sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device', screenPrimaryChannelObject.screenPrimarySpd);
SetSensorColorSpace(screenCalObj, colorDirectionParams.T_receptors, colorDirectionParams.S);

%% Create ISETBio display from the calibration file.
[ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj,'verbose',VERBOSE);

%% Set up the background screen primaries.
backgroundScreenPrimaryObject = SetupBackground(colorDirectionParams,screenCalObj,backgroundChannelObject,'verbose',VERBOSE);

%% Make a monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image. Multiply by the max
% contrast vector to get the LMS contrast image. The function includes the
% quantization of the gabor image.
[rawMonochromeUnquantizedContrastGaborImage, rawMonochromeUnquantizedContrastGaborCal, rawMonochromeContrastGaborCal, ...
    stimulusN, centerN, stimulusHorizSizeDeg, stimulusHorizSizeMeters] = ...
      MakeMonochromeContrastGaborMel(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,'verbose',VERBOSE,'nQuantizeBits',nQuantizeBits);

%% Get cone contrast/excitation gabor image.
%
% The SpectralCalISETBio routine encapsulated the inline code below more, but the called
% functions did enough else that I didn't want to merge it right in at this juncture.  I did
% set up the encapsulating structures at the top level to match, including making them
% cell arrays.  But here the sf and contrast variables are fixed for now, and just return
% one entry. To return to.
%
% Compute desired contrast
ss = 1; cc = 1;
standardGaborCalObject.desiredContrastGaborCal{ss,cc} = colorDirectionParams.spatialGaborTargetContrast*colorDirectionParams.targetStimulusContrastDir*rawMonochromeContrastGaborCal{ss};

% Convert cone contrast to excitations
standardGaborCalObject.desiredExcitationsGaborCal{ss,cc} = ContrastToExcitation(standardGaborCalObject.desiredContrastGaborCal{ss,cc},backgroundScreenPrimaryObject.screenBgExcitations);

% Get primaries using standard calibration code, and desired spd without
% quantizing.
standardGaborCalObject.standardPrimariesGaborCal{ss,cc} = SensorToPrimary(screenCalObj,standardGaborCalObject.desiredExcitationsGaborCal{ss,cc});
standardGaborCalObject.desiredSpdGaborCal{ss,cc} = PrimaryToSpd(screenCalObj,standardGaborCalObject.standardPrimariesGaborCal{ss,cc});

% Gamma correct and quantize (if gamma method set to 2 above; with gamma
% method set to zero there is no quantization).  Then convert back from
% the gamma corrected settings.
standardGaborCalObject.standardSettingsGaborCal{ss,cc} = PrimaryToSettings(screenCalObj,standardGaborCalObject.standardPrimariesGaborCal{ss,cc});
standardGaborCalObject.standardPredictedPrimariesGaborCal{ss,cc} = SettingsToPrimary(screenCalObj,standardGaborCalObject.standardSettingsGaborCal{ss,cc});
standardGaborCalObject.standardPredictedExcitationsGaborCal{ss,cc} = PrimaryToSensor(screenCalObj,standardGaborCalObject.standardPredictedPrimariesGaborCal{ss,cc});
standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc} = ExcitationsToContrast(standardGaborCalObject.standardPredictedExcitationsGaborCal{ss,cc},backgroundScreenPrimaryObject.screenBgExcitations);

% Plot of how well standard method does in obtaining desired contrast
if (VERBOSE)
    standardContrastScatterFig = figure;
    set(gcf,'Position',[100 100 1500 600]);
    subplot(1,4,1);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(1,:),standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(1,:),'r+');
    fprintf('Standard image max L contrast: %0.3f\n',max(abs(standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(1,:))));
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired L contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted L contrast','FontName','Helvetica','FontSize',18);
    title('Standard image method','FontName','Helvetica','FontSize',18);

    subplot(1,4,2);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(2,:),standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(2,:),'g+');
    fprintf('Standard image max M contrast: %0.3f\n',max(abs(standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(2,:))));
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired M contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted M contrast','FontName','Helvetica','FontSize',18);
    title('Standard image method','FontName','Helvetica','FontSize',18);

    subplot(1,4,3);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(3,:),standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(3,:),'b+');
    fprintf('Standard image max S contrast: %0.3f\n',max(abs(standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(3,:))));
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired S contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted S contrast','FontName','Helvetica','FontSize',18);
    title('Standard image method','FontName','Helvetica','FontSize',18);

    subplot(1,4,4);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(4,:),standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(4,:),'c+');
    fprintf('Standard image max Mel contrast: %0.3f\n',max(abs(standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc}(4,:))));
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired MEL contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted MEL contrast','FontName','Helvetica','FontSize',18);
    title('Standard image method','FontName','Helvetica','FontSize',18);

    saveas(standardContrastScatterFig,fullfile(outputDir,'standardContrastScatterFig.tiff'),'tiff');
end

%% Set up table of contrasts for all possible settings
[ptCldSettingsCal, ptCldContrastCal] = SetupContrastPointLookup(screenCalObj,backgroundScreenPrimaryObject.screenBgExcitations,'verbose',VERBOSE);

%% Get image from table, in cal format
uniqueQuantizedSettingsGaborCal = SettingsFromLookup(standardGaborCalObject.desiredContrastGaborCal{ss,cc},ptCldContrastCal,ptCldSettingsCal);

% Print out min/max of settings
fprintf('Gabor image min/max settings: %0.3f, %0.3f\n',min(uniqueQuantizedSettingsGaborCal(:)), max(uniqueQuantizedSettingsGaborCal(:)));

% Get contrasts we think we have obtianed
uniqueQuantizedExcitationsGaborCal = SettingsToSensor(screenCalObj,uniqueQuantizedSettingsGaborCal);
uniqueQuantizedContrastGaborCal = ExcitationsToContrast(uniqueQuantizedExcitationsGaborCal,backgroundScreenPrimaryObject.screenBgExcitations);

%% Convert representations we want to take forward to image format.
%
% Note that we set the settings field to the uniqueQuantized settings,
% because these are what we want to use.  This differs from what
% is done in the routine MakeImageSettingsFromPtCld.  I think that
% routine is in error, but not changing pending response from Semin.
gaborImageObject.desiredContrastGaborImage{ss,cc} = CalFormatToImage(standardGaborCalObject.desiredContrastGaborCal{ss,cc},stimulusN,stimulusN);
gaborImageObject.settingsGaborImage{ss,cc} = CalFormatToImage(uniqueQuantizedSettingsGaborCal,stimulusN,stimulusN);
gaborImageObject.excitationsGaborImage{ss,cc} = CalFormatToImage(uniqueQuantizedExcitationsGaborCal,stimulusN,stimulusN);
gaborImageObject.contrastGaborImage{ss,cc} = CalFormatToImage(uniqueQuantizedContrastGaborCal,stimulusN,stimulusN);
lContrastNominal = max(abs(uniqueQuantizedContrastGaborCal(1,:)));
mContrastNominal = max(abs(uniqueQuantizedContrastGaborCal(2,:)));
sContrastNominal = max(abs(uniqueQuantizedContrastGaborCal(3,:)));
melContrastNominal = max(abs(uniqueQuantizedContrastGaborCal(4,:)));
fprintf('Quantized unique point image max L contrast: %0.3f\n',lContrastNominal );
fprintf('Quantized unique point image max M contrast: %0.3f\n',mContrastNominal);
fprintf('Quantized unique point image max S contrast: %0.3f\n',sContrastNominal);
fprintf('Quantized unique point image max Mel contrast: %0.3f\n',melContrastNominal);
lmsmelContrastNominal = [lContrastNominal ; mContrastNominal ; sContrastNominal; melContrastNominal];

% Plot of how well point cloud method does in obtaining desired contratsfigure; clf;
if (VERBOSE)
    pointCloutContrastScatterFig = figure;
    set(gcf,'Position',[100 100 1500 600]);
    subplot(1,4,1);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(1,:),uniqueQuantizedContrastGaborCal(1,:),'r+');
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired L contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted L contrast','FontName','Helvetica','FontSize',18);
    title('Nominal vs Desired','FontName','Helvetica','FontSize',18);

    subplot(1,4,2);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(2,:),uniqueQuantizedContrastGaborCal(2,:),'g+');
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired M contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted M contrast','FontName','Helvetica','FontSize',18);
    title('Nominal vs Desired','FontName','Helvetica','FontSize',18);

    subplot(1,4,3);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(3,:),uniqueQuantizedContrastGaborCal(3,:),'b+');
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired S contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted S contrast','FontName','Helvetica','FontSize',18);
    title('Nominal vs Desired','FontName','Helvetica','FontSize',18);

    subplot(1,4,4);
    set(gca,'FontName','Helvetiaca','FontSize',16);
    plot(standardGaborCalObject.desiredContrastGaborCal{ss,cc}(4,:),uniqueQuantizedContrastGaborCal(4,:),'c+');
    axis('square');
    xlim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]); ylim([-1.2*spatialGaborTargetContrast 1.2*spatialGaborTargetContrast]);
    xlabel('Desired MEL contrast','FontName','Helvetica','FontSize',18);
    ylabel('Predicted MEL contrast','FontName','Helvetica','FontSize',18);
    title('Nominal vs Desired','FontName','Helvetica','FontSize',18);

    saveas(pointCloutContrastScatterFig,fullfile(outputDir,'lookupContrastScatterFig.tiff'),'tiff');
end

%% Put the image into an ISETBio scene.
% 
% As I read the subroutine, it is using the standard rather than lookup
% table data here. That's fine, but we should stay aware of it.
ISETBioGaborObject = MakeISETBioSceneFromImageMel(colorDirectionParams,gaborImageObject, ...
    ISETBioDisplayObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg,'verbose',false);

% Go back to the RGB image starting with the ISETBio representation.
%
% This check commented out because it is written assuming that the ISETBio
% scene settings were computed using the standard method, and we are using
% the exhaustive search settings for improved precision.  Leaving it here
% to note that we could rewrite the check if we wanted to, but it would be
% slow.
%
% fromISETBioGaborCalObject = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborObject,standardGaborCalObject,'verbose',VERBOSE);
% t1 = ImageToCalFormat(gaborImageObject.excitationsGaborImage{ss,cc});
% t2 = fromISETBioGaborCalObject.settingsFromISETBioGaborCal;
% diff = max(abs(t1(:)-t2(:)));
% if (diff > 1e-6)
%     error('ISETBio scene does not reproduce settings\n');
% end

%% SRGB image via XYZ, scaled to display
predictedXYZCal = colorDirectionParams.T_xyz * ISETBioGaborObject.ISETBioGaborImageSpdCal{ss,cc};
SRGBPrimaryCal = XYZToSRGBPrimary(predictedXYZCal);
scaleFactor = max(SRGBPrimaryCal(:));
SRGBCal = SRGBGammaCorrect(SRGBPrimaryCal/(2*scaleFactor),0);
SRGBImage = uint8(CalFormatToImage(SRGBCal,stimulusN,stimulusN));

% Show the SRGB image
figure; imshow(SRGBImage);
title('SRGB Gabor Image');
imwrite(SRGBImage,fullfile(outputDir,'sRGBImage.tiff'),'tiff');

%% L plane image as grayscale sRGB
predictedLPlaneCal = colorDirectionParams.T_receptors(1,:) * ISETBioGaborObject.ISETBioGaborImageSpdCal{ss,cc};
predictedLPlaneSRGBPlaneCal = [predictedLPlaneCal ; predictedLPlaneCal ; predictedLPlaneCal];
SRGBLPlanePrimaryCal = predictedLPlaneSRGBPlaneCal;
scaleFactor = mean(predictedLPlaneSRGBPlaneCal(:));
SRGBLPlaneCal = SRGBGammaCorrect(SRGBLPlanePrimaryCal/(2*scaleFactor),0);
SRGBLPlaneImage = uint8(CalFormatToImage(SRGBLPlaneCal,stimulusN,stimulusN));
figure; imshow(SRGBLPlaneImage);
title('SRGB L Plane Gabor Image');
imwrite(SRGBLPlaneImage,fullfile(outputDir,'sRGBLPlaneImage.tiff'),'tiff');

%% M plane image as grayscale sRGB
predictedMPlaneCal = colorDirectionParams.T_receptors(2,:) * ISETBioGaborObject.ISETBioGaborImageSpdCal{ss,cc};
predictedMPlaneSRGBPlaneCal = [predictedMPlaneCal ; predictedMPlaneCal ; predictedMPlaneCal];
SRGBMPlanePrimaryCal = predictedMPlaneSRGBPlaneCal;
scaleFactor = mean(predictedMPlaneSRGBPlaneCal(:));
SRGBMPlaneCal = SRGBGammaCorrect(SRGBMPlanePrimaryCal/(2*scaleFactor),0);
SRGBMPlaneImage = uint8(CalFormatToImage(SRGBMPlaneCal,stimulusN,stimulusN));
figure; imshow(SRGBMPlaneImage);
title('SRGB M Plane Gabor Image');
imwrite(SRGBMPlaneImage,fullfile(outputDir,'sRGBMPlaneImage.tiff'),'tiff');

%% S plane image as grayscale sRGB
predictedSPlaneCal = colorDirectionParams.T_receptors(3,:) * ISETBioGaborObject.ISETBioGaborImageSpdCal{ss,cc};
predictedSPlaneSRGBPlaneCal = [predictedSPlaneCal ; predictedSPlaneCal ; predictedSPlaneCal];
SRGBSPlanePrimaryCal = predictedSPlaneSRGBPlaneCal;
scaleFactor = mean(predictedSPlaneSRGBPlaneCal(:));
SRGBSPlaneCal = SRGBGammaCorrect(SRGBSPlanePrimaryCal/(2*scaleFactor),0);
SRGBSPlaneImage = uint8(CalFormatToImage(SRGBSPlaneCal,stimulusN,stimulusN));
figure; imshow(SRGBSPlaneImage);
title('SRGB S Plane Gabor Image');
imwrite(SRGBSPlaneImage,fullfile(outputDir,'sRGBSPlaneImage.tiff'),'tiff');

%% Mel plane image as grayscale sRGB
predictedMelPlaneCal = colorDirectionParams.T_receptors(4,:) * ISETBioGaborObject.ISETBioGaborImageSpdCal{ss,cc};
predictedMelPlaneSRGBPlaneCal = [predictedMelPlaneCal ; predictedMelPlaneCal ; predictedMelPlaneCal];
SRGBMelPlanePrimaryCal = predictedMelPlaneSRGBPlaneCal;
scaleFactor = mean(predictedMelPlaneSRGBPlaneCal(:));
SRGBMelPlaneCal = SRGBGammaCorrect(SRGBMelPlanePrimaryCal/(2*scaleFactor),0);
SRGBMelPlaneImage = uint8(CalFormatToImage(SRGBMelPlaneCal,stimulusN,stimulusN));
figure; imshow(SRGBMelPlaneImage);
title('SRGB Mel Plane Gabor Image');
imwrite(SRGBMelPlaneImage,fullfile(outputDir,'sRGBMelPlaneImage.tiff'),'tiff');

%% Show the settings image
figure; clf;
imshow(gaborImageObject.settingsGaborImage{ss,cc});
title('Image of settings');
imwrite(gaborImageObject.settingsGaborImage{ss,cc},fullfile(outputDir,'settingsImage.tiff'),'tiff');

%% Plot slice through predicted LMS and Mel contrast image.
%
% Set the plot limit axis.
stimulusN = size(gaborImageObject.contrastGaborImage{ss,cc},1);
stimulusDegs = (stimulusSizeDeg/2)*((1:stimulusN)-stimulusN/2)/(stimulusN/2);
plotAxisLimit = 1.2*100 * colorDirectionParams.spatialGaborTargetContrast;
PlotSliceContrastGaborImage(gaborImageObject.contrastGaborImage{ss,cc}, gaborImageObject.contrastGaborImage{ss,cc}, ...
    'plotAxisLimit', plotAxisLimit, 'verbose', VERBOSE, 'xAxisValues', stimulusDegs, 'lineWidth',2,'markerSize',0);
set(gca,'FontName','Helvetiaca','FontSize',16);
%title('LMS and Mel Contrast','FontName','Helvetiaca','FontSize',14);
legend({sprintf('L: %0.1f%%',100*lContrastNominal), sprintf('M: %0.1f%%',100*mContrastNominal), ...
    sprintf('S: %0.1f%%',100*sContrastNominal), sprintf('Mel: %0.1f%%',100*melContrastNominal)},'FontName','Helvetiaca','FontSize',12,'Location','SouthEast');
xlabel('Relative Position (deg)','FontName','Helvetiaca','FontSize',18);
ylabel('Nominal Contrast','FontName','Helvetiaca','FontSize',18);
saveas(gcf,fullfile(outputDir,'LMSMelGaborSlice.tiff'),'tiff');
% tempCal = ImageToCalFormat(gaborImageObject.contrastGaborImage{ss,cc});
% max(tempCal,[],2)

%% Plot the settings
plotAxisLimit = 255;
PlotSliceContrastGaborImage(gaborImageObject.settingsGaborImage{ss,cc}, gaborImageObject.settingsGaborImage{ss,cc}, ...
    'plotAxisLimit', plotAxisLimit, 'plotAxisPos', true, 'verbose', VERBOSE, 'xAxisValues', stimulusDegs, 'lineWidth',2,'markerSize',0);
set(gca,'FontName','Helvetiaca','FontSize',16);
%title('Primary Settings','FontName','Helvetiaca','FontSize',14);
legend({'Primary 1', 'Primary 2', 'Primary 3'},'FontName','Helvetiaca','FontSize',12,'Location','NorthEast');
xlabel('Relative Position (deg)','FontName','Helvetiaca','FontSize',18);
ylabel('Primary Settings','FontName','Helvetiaca','FontSize',18);
saveas(gcf,fullfile(outputDir,'PrimarySettingsGaborSlice.tiff'),'tiff');

%% Generate some settings values corresponding to known contrasts
%
% The reason for this is to measure and check these.  This logic follows
% how we handled an actual gabor image above.
rawMonochromeUnquantizedContrastCheckCal = [0 0.05 -0.05 0.10 -0.10 0.15 -0.15 0.20 -0.20 0.25 -0.25 0.5 -0.5 1 -1];
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
desiredContrastCheckCal = colorDirectionParams.spatialGaborTargetContrast*colorDirectionParams.targetStimulusContrastDir*rawMonochromeContrastCheckCal;
desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,backgroundScreenPrimaryObject.screenBgExcitations);

% For each check calibration find the settings that
% come as close as possible to producing the desired excitations.
%
% If we measure for a uniform field the spectra corresopnding to each of
% the settings in the columns of ptCldScreenSettingsCheckCall, then
% compute the cone contrasts with respect to the backgound (0 contrast
% measurement, first settings), we should approximate the cone contrasts in
% desiredContrastCheckCal. 
ptCldScreenSettingsCheckCal = SettingsFromLookup(desiredContrastCheckCal,ptCldContrastCal,ptCldSettingsCal);
ptCldScreenPrimariesCheckCal = SettingsToPrimary(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenSpdCheckCal = PrimaryToSpd(screenCalObj,ptCldScreenPrimariesCheckCal);
ptCldScreenExcitationsCheckCal = SettingsToSensor(screenCalObj,ptCldScreenSettingsCheckCal);
ptCldScreenContrastCheckCal = ExcitationsToContrast(ptCldScreenExcitationsCheckCal,backgroundScreenPrimaryObject.screenBgExcitations);
figure; clf; hold on;
plot(desiredContrastCheckCal(4,:),ptCldScreenContrastCheckCal(4,:),'co','MarkerSize',10,'MarkerFaceColor','c');
plot(desiredContrastCheckCal(3,:),ptCldScreenContrastCheckCal(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');
plot(desiredContrastCheckCal(2,:),ptCldScreenContrastCheckCal(2,:),'go','MarkerSize',10,'MarkerFaceColor','g');
plot(desiredContrastCheckCal(1,:),ptCldScreenContrastCheckCal(1,:),'ro','MarkerSize',10,'MarkerFaceColor','r');
xlim([0 plotAxisLimit/100]); ylim([0 plotAxisLimit/100]); axis('square');
xlabel('Desired'); ylabel('Obtained');
title('Desired versus obtained check contrasts');
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
    screenPrimarySpdCheck(:,pp) = PrimaryToSpd(channelCalObjs{pp},SettingsToPrimary(channelCalObjs{pp}, screenPrimaryChannelObject.screenPrimarySettings(:,pp)));
end
figure; clf; hold on
set(gca,'FontName','Helvetiaca','FontSize',16);
plot(colorDirectionParams.wls, screenPrimaryChannelObject.screenPrimarySpd(:,1),'r','LineWidth',4);
plot(colorDirectionParams.wls, screenPrimaryChannelObject.screenPrimarySpd(:,2),'g','LineWidth',4);
plot(colorDirectionParams.wls, screenPrimaryChannelObject.screenPrimarySpd(:,3),'b','LineWidth',4);
plot(colorDirectionParams.wls, screenPrimarySpdCheck,'k','LineWidth',2);
xlabel('Wavelength','FontName','Helvetiaca','FontSize',18);
ylabel('Radiance','FontName','Helvetiaca','FontSize',18);
%title('Primary spds','FontName','Helvetiaca','FontSize',20);
legend({'Primary 1', 'Primary 2', 'Primary 3'},'FontName','Helvetiaca','FontSize',12,'Location','NorthWest');
saveas(gcf,fullfile(outputDir,'primarySpds.tiff'),'tiff');

%% Save out what we need 
save(sceneOutputFilename,'lmsmelContrastNominal', 'colorDirectionParams','screenPrimaryChannelObject','backgroundChannelObject','backgroundScreenPrimaryObject', ...
    'ISETBioDisplayObject','screenSizeObject','screenCalObjFromISETBio', ...
    'ISETBioGaborObject', ...
    '-v7.3');
disp('Data has been saved successfully!');
