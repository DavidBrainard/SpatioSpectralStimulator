% SpectralTestCal
%
% Start exploring spectral fits with swubprimarys, this
% version using the calibration structures.
%
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Define calibration filenames/params.
%
% This is a standard calibration file for the DLP projector,
% with the subprimaries set to something.  As we'll see below,
% we're going to rewrite those.
projectorCalName = 'SACC';
projectorNInputLevels = 256;

% These are the calibration files for each of the primaries, which
% then entails measuring the spectra of all the subprimaries for that
% primary.
subprimaryCalNames = {'SACCPrimary1' 'SACCPrimary1' 'SACCPrimary1'};
subprimaryNInputLevels = 253;

%% Load projector calibration.
projectorCal = LoadCalFile(projectorCalName);
projectorCalObj = ObjectToHandleCalOrCalStruct(projectorCal);
CalibrateFitGamma(projectorCalObj, projectorNInputLevels);

%% Load subprimary calibrations.
nPrimaries = 3;
subprimaryCals = cell(nPrimaries ,1);
subprimaryCalObjs = cell(nPrimaries ,1);
for cc = 1:length(subprimaryCalNames)
    subprimaryCals{cc} = LoadCalFile(subprimaryCalNames{cc});

    subprimaryCalObjs{cc} = ObjectToHandleCalOrCalStruct(subprimaryCals{cc});
    CalibrateFitGamma(subprimaryCalObjs{cc}, subprimaryNInputLevels);
end

%% Get out some data to work with.
%
% This is from the subprimary calibration file.
S = subprimaryCalObjs{1}.get('S');
wls = SToWls(S);
ambientSpd = subprimaryCalObjs{1}.get('P_ambient');
if (isempty(ambientSpd))
    subprimaryCalObjs{1}.P_ambient = zeros(size(wls));
end
P_device = subprimaryCalObjs{1}.get('P_device');
gammaInput = subprimaryCalObjs{1}.get('gammaInput');
gammaTable = subprimaryCalObjs{1}.get('gammaTable');
gammaMeasurements = subprimaryCals{1}.rawData.gammaCurveMeanMeasurements;
[nSubprimaries,nMeas,~] = size(gammaMeasurements);

%% Cone fundamentals and XYZ CMFs.
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
T_cones = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);
load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,S);

%% Let's look at little at the subprimary calibration.
%
% Eventually this will be handled by the analyze program,
% when it is generalized for more than three primaries.  But
% we are impatient people so we will hack something up here.
PLOT_SUBPRIMARYINVARIANCE = false;
if (PLOT_SUBPRIMARYINVARIANCE)
    for pp = 1:nSubprimaries
        maxSpd = squeeze(gammaMeasurements(pp,end,:));
        figure;
        subplot(1,2,1); hold on;
        plot(wls,maxSpd,'r','LineWidth',3);
        for mm = 1:nMeas-1
            temp = squeeze(gammaMeasurements(pp,mm,:));
            plot(wls,temp,'k','LineWidth',1);
        end
        subplot(1,2,2); hold on
        plot(wls,maxSpd,'r','LineWidth',3);
        for mm = 1:nMeas-1
            temp = squeeze(gammaMeasurements(pp,mm,:));
            scaleFactor = temp\maxSpd;
            plot(wls,scaleFactor*temp,'k','LineWidth',1);
        end
    end
end

%% Plot subprimary gamma functions.
PLOT_SUBPRIMARYGAMMA = false;
if (PLOT_SUBPRIMARYGAMMA)
    for pp = 1:nSubprimaries
        figure; hold on;
        plot(subprimaryCals{1}.rawData.gammaInput,subprimaryCals{1}.rawData.gammaTable(:,pp),'ko','MarkerSize',12,'MarkerFaceColor','k');
        plot(gammaInput,gammaTable(:,pp),'k','LineWidth',2);
    end 
end

%% Plot x,y if desired.
PLOT_SUBPRIMARYCHROMATICITY = false;
if (PLOT_SUBPRIMARYCHROMATICITY)
    figure; hold on;
    for pp = 1:nSubprimaries
        for mm = 1:nMeas
            % XYZ calculation for each measurement
            spd_temp = squeeze(gammaMeasurements(pp,mm,:));      
            XYZ_temp = T_xyz*spd_temp; 
            xyY_temp = XYZToxyY(XYZ_temp);
            
            plot(xyY_temp(1,:),xyY_temp(2,:),'r.','Markersize',10); % Coordinates of the subprimary
            xlabel('CIE x');
            ylabel('CIE y');
        end
    end
    
    % Add spectrum locus to plot, connected end to end
    colorgamut=XYZToxyY(T_xyz); 
    colorgamut(:,end+1)=colorgamut(:,1);
    plot(colorgamut(1,:),colorgamut(2,:),'k-'); 
end

%% Background xy.
%
% Specify the chromaticity, but we'll chose the luminance based
% on the range available in the device.
targetBgxy = [0.3127 0.3290]';

%% Target color direction and max contrasts.
%
% This is the basic desired modulation direction positive excursion. We go
% equally in positive and negative directions.
targetLMSContrast = [1 -1 0]';

%% Specify desired primary properties.
%
% These are the target contrasts for the three primaries. We want these to
% span a triangle around the line specified above. Here we define that
% triangle by hand.  May need a little fussing for other directions, and
% might be able to autocompute good choices.
target1MaxLMSContrast = [-1 1 0]';
target2MaxLMSContrast = [1 -1 0.5]';
target3MaxLMSContrast = [1 -1 -0.5]';

% We may not need the whole direction excursion above. The first number is
% the amount we want to use, the second has a little headroom so we don't
% run into numerical error at the edges. The second number is used when
% defining the three primaries, the first when computing desired weights on
% the primaries.
targetContrastReMax = 0.05;
targetPrimaryHeadroom = 1.1;
targetContrastReMaxWithHeadroom = targetPrimaryHeadroom*targetContrastReMax;
plotAxisLimit = 2;

%% Comment this better later on.
%
% When we compute a specific image, we may not want full contrast available
% with the primaries. This tells us fraction of max available relative to
% ledContrastReMax.
imageModulationContrast = 0.05/targetContrastReMax;

%% Image spatial parameters.
sineFreqCyclesPerImage = 6;
gaborSdImageFraction = 0.1;

% Image size in pixels
imageN = 512;

%% Computational bit depth.
%
% This is a computational bit depth that we use to define the lookup table
% between contrast and primary values.
fineBits = 14;
nFineLevels = 2^fineBits;

%% Get half on spectrum.
%
% This is useful for scaling things reasonably - we start with half of the
% available range of the primaries.
halfOnSubprimaries = 0.5*ones(nSubprimaries,1);
halfOnSpd = PrimaryToSpd(subprimaryCalObjs{1},halfOnSubprimaries);

%% Make sure gamma correction behaves well with unquantized conversion.
SetGammaMethod(subprimaryCalObjs{1},0);
halfOnSettings = PrimaryToSettings(subprimaryCalObjs{1},halfOnSubprimaries);
halfOnPrimariesChk = SettingsToPrimary(subprimaryCalObjs{1},halfOnSettings);
if (max(abs(halfOnSubprimaries-halfOnPrimariesChk)) > 1e-8)
    error('Gamma self-inversion not sufficiently precise');
end

%% Use quantized conversion from here on.
%
% Comment in the line that refits the gamma to see
% effects of extreme quantization below
%
% CalibrateFitGamma(subprimaryCalObjs{1},10);
SetGammaMethod(subprimaryCalObjs{1},2);
SetGammaMethod(subprimaryCalObjs{2},2);
SetGammaMethod(subprimaryCalObjs{3},2);

%% Use extant machinery to get primaries from spectrum.
%
% This isn't used in our calculations.  Any difference in the
% two lines here reflects a bug in the SpdToPrimary/PrimaryToSpd pair.  
halfOnPrimariesChk = SpdToPrimary(subprimaryCalObjs{1},halfOnSpd);
halfOnSpdChk = PrimaryToSpd(subprimaryCalObjs{1},halfOnPrimariesChk);
figure; hold on;
plot(wls,halfOnSpd,'r','LineWidth',3);
plot(wls,halfOnSpdChk,'k','LineWidth',1);

%% Show effect of quantization.
%
% It's very small at the nominal 252 levels of the subprimaries, but will
% increase if you refit the gamma functios to a small number of levels.
halfOnPrimariesChk = SpdToPrimary(subprimaryCalObjs{1},halfOnSpd);
halfOnSettingsChk = PrimaryToSettings(subprimaryCalObjs{1},halfOnPrimariesChk);
halfOnPrimariesChk1 = SettingsToPrimary(subprimaryCalObjs{1},halfOnSettingsChk);
halfOnSpdChk1 = PrimaryToSpd(subprimaryCalObjs{1},halfOnPrimariesChk1);
plot(wls,halfOnSpdChk1,'g','LineWidth',1);

%% Set up basis to try to keep spectra close to.
%
% This is how we enforce a smoothness or other constraint
% on the spectra.
basisType = 'fourier';
nFourierBases = 7;
switch (basisType)
    case 'cieday'
        load B_cieday
        B_natural = SplineSpd(S_cieday,B_cieday,S);
    case 'fourier'
        B_natural = MakeFourierBasis(S,nFourierBases);
    otherwise
        error('Unknown basis set specified');
end

% Define wavelength range that will be used to enforce the smoothnes
% thorugh the projection onto an underlying basis set.  We don't the whole
% visible spectrum as putting weights on the extrema where people are not
% sensitive costs us smoothness in the spectral region we care most about.
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(wls > lowProjectWl & wls < highProjectWl);

%% Find background primaries to acheive desired xy at intensity scale of display.
% Set parameters for getting desired background primaries.
primaryHeadRoom = 0;
targetLambda = 3;
targetBgXYZ = xyYToXYZ([targetBgxy ; 1]);

% Make a loop for getting background primaries for all primaries.
for pp = 1:nPrimaries
    [bgPrimaries(:,pp),obtainedBgSpd(:,pp),obtainedBgXYZ(:,pp)] = FindDesiredBackgroundPrimaries(targetBgXYZ,T_xyz,subprimaryCalObjs{pp}, ...
        B_natural,projectIndices,primaryHeadRoom,targetLambda,'Scale',true,'Verbose',true);
end

if (any(bgPrimaries < 0) | any(bgPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end

%% Find target primaries with desired LMS contrast.
% Set parameters for getting desired target primaries.
targetMaxLMSContrast = [target1MaxLMSContrast target2MaxLMSContrast target3MaxLMSContrast];
targetContrastReMax = 0.05;
targetPrimaryHeadroom = 1.1;
primaryHeadroom = 0;
targetLambda = 3;

% Make a loop for getting isolating primaries for all primaries.
for pp = 1:nPrimaries
    otherPrimaries = setdiff(1:nPrimaries,pp);
    extraAmbientSpd = 0;
    % Set extra ambient spd.
    for oo = 1:length(otherPrimaries)
        extraAmbientSpd = extraAmbientSpd + obtainedBgSpd(:,otherPrimaries(oo));
    end
    % Get isolating primaries.
    [isolatingPrimaries(:,pp),isolatingPrimariesQuantized(:,pp),isolatingSpd(pp,:),isolatingContrast(pp,:)] = FindDesiredContrastTargetPrimaries(targetMaxLMSContrast(:,pp), ...
        targetPrimaryHeadroom,targetContrastReMax,bgPrimaries(:,pp), ...
        T_cones,subprimaryCalObjs{pp},B_natural,projectIndices,primaryHeadroom,targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
end


%% Measure the desired target primaries (THIS PART HAS BEEN UPDATED - SEMIN)
% This result will be used to compute the image.
%
% Make a loop for measuring all primaries.
for pp = 1:nPrimaries
    [isolatingSpdMeasured(pp,:)] = MeasureDesiredTargetPrimaries(isolatingPrimaries(:,pp),subprimaryNInputLevels,subprimaryCalObjs{pp},pp,'projectorMode',true,'measurementOption',false,'verbose',false);
end

%% How close are spectra to subspace defined by basis?
theBgNaturalApproxSpd = B_natural*(B_natural(projectIndices,:)\obtainedBgSpd(projectIndices));
isolatingNaturalApproxSpd1 = B_natural*(B_natural(projectIndices,:)\isolatingSpd(1,projectIndices)');
isolatingNaturalApproxSpd2 = B_natural*(B_natural(projectIndices,:)\isolatingSpd(2,projectIndices)');
isolatingNaturalApproxSpd3 = B_natural*(B_natural(projectIndices,:)\isolatingSpd(3,projectIndices)');

% Plot
figure; clf;
subplot(2,2,1); hold on
plot(wls,obtainedBgSpd,'b','LineWidth',2);
plot(wls,theBgNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),obtainedBgSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theBgNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Background');
%ylim([0 2]);
subplot(2,2,2); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(1,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd1,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(1,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd1(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 1');
%ylim([0 2]);
subplot(2,2,3); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(2,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd2,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(2,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd2(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 2');
%ylim([0 2]);
subplot(2,2,4); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(3,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd3,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(3,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd3(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 3');
%ylim([0 2]);

%% Set the projector primaries
%
% We want these to match those we set up with the
% subprimary calculations above.  Need to reset
% sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
projectorCalObj.set('P_device',isolatingSpd');
SetSensorColorSpace(projectorCalObj,T_cones,S);

%% Set gamma method
%
% If we set to 0, there is no quantization and the result is excellent.
% If we set to 2, this is quantized at 256 levels and the result is a mess.
% This is what the nearest neighbor code, which I have stripped out, fixes.
% Need to put that back in, but maybe not in the same way.
SetGammaMethod(projectorCalObj,0);

%% Make Gabor patch in range 0-1.
%
% This is our contrast modulation image.  Multiply
% by the max contrast vector to get the LMS contrast image.
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawMonochromeSineImage = MakeSineImage(0,sineFreqCyclesPerImage,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawMonochromeGaborImage = imageModulationContrast*rawMonochromeSineImage.*gaussianWindow;
rawMonochromeGaborCal = ImageToCalFormat(rawMonochromeGaborImage);

fprintf('Making Gabor LMS contrast image\n');
bgSpd = sum(obtainedBgSpd,2);
bgLMS = T_cones*bgSpd;
theDesiredLMSContrastGaborCal = targetContrastReMax*targetLMSContrast*rawMonochromeGaborCal;
theDesiredLMSExcitationGaborCal = ContrastToExcitation(theDesiredLMSContrastGaborCal,bgLMS);
thePrimariesGaborCal = SensorToPrimary(projectorCalObj,theDesiredLMSExcitationGaborCal);
theSettingsGaborCal = PrimaryToSettings(projectorCalObj,thePrimariesGaborCal);
thePredictedPrimariesGaborCal = SettingsToPrimary(projectorCalObj,theSettingsGaborCal);
thePredictedLMSExcitationsGaborCal = PrimaryToSensor(projectorCalObj,thePredictedPrimariesGaborCal);
thePredictedSpdGaborCal = PrimaryToSpd(projectorCalObj,thePredictedPrimariesGaborCal);

meanLMS = mean(thePredictedLMSExcitationsGaborCal,2);
theDesiredLMSContrastGaborImage = CalFormatToImage(theDesiredLMSContrastGaborCal,imageN,imageN);
thePredictedLMSContrastCal = ExcitationsToContrast(thePredictedLMSExcitationsGaborCal,meanLMS);
thePredictedLMSContrastImage = CalFormatToImage(thePredictedLMSContrastCal,imageN,imageN);
theSettingsGaborImage = CalFormatToImage(theSettingsGaborCal,imageN,imageN);

%% Set up exhaustive point cloud for all possible quantized input settings
allSettingsCal = zeros(3,256^3);
idx = 1;
for ii = 0:255
    for jj = 0:255
        for kk = 0:255
            allSettingsCal(:,idx) = [ii jj kk]';
            idx = idx+1;
        end
    end
end
allSettingsCal = allSettingsCal/255;
allSensorCal = SettingsToSensor(projectorCalObj,allSettingsCal);
allSensorContrastCal = ExcitationsToContrast(thePointCloudLMSExcitationsGaborCal,meanLMS);
allSensorPtCloud = pointCloud(allSensorContrastCal');

%% Do SensorToSettings by finding nearest points in point cloud
thePointCloudSettingsGaborCal = zeros(3,size(theDesiredLMSContrastGaborCal,2));
for ll = 1:size(theDesiredLMSContrastGaborCal,2) 
    minIndex = findNearestNeighbors(allSensorPtCloud,theDesiredLMSContrastGaborCal(:,ll)',1);
    thePointCloudSettingsGaborCal(:,ll) = allSensorCal(:,minIndex);
end
thePointCloudLMSExcitationsGaborCal = SettingsToSensor(projectorCalObj,thePointCloudSettingsGaborCal);
thePointCloudLMSContrastGaborCal = ExcitationsToContrast(thePointCloudLMSExcitationsGaborCal,meanLMS);
thePointCloudLMSContrastGaborImage = CalFormatToImage(thePointCloudLMSContrastGaborCal,imageN,imageN);

% SRGB image via XYZ
thePredictedXYZCal = T_xyz*thePredictedSpdGaborCal;
theSRGBPrimaryCal = XYZToSRGBPrimary(thePredictedXYZCal);
scaleFactor = max(theSRGBPrimaryCal(:));
theSRGBCal = SRGBGammaCorrect(theSRGBPrimaryCal/(2*scaleFactor),0);
theSRGBImage = uint8(CalFormatToImage(theSRGBCal,imageN,imageN));

% Show the SRGB image
figure; imshow(theSRGBImage);

% Show the settings image
figure; clf;
imshow(theSettingsGaborImage);

%% Show this image on the DLP, and it should look more or less like
% the sRGB image we display below.
projectorSettingsImage = theSettingsGaborImage;
testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
testFilename = fullfile(testFiledir,'testImageData1');
save(testFilename,'projectorSettingsImage','isolatingPrimaries');

%% Plot slice through predicted LMS contrast image.
figure; hold on
plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Plot slice through point cloud LMS contrast image.
figure; hold on
plot(1:imageN,100*thePointCloudLMSContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*thePointCloudLMSContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*thePointCloudLMSContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Light level tests.
%
% PupilDiameter
pupilDiameterMM = 4;
theStimulusExtentDeg = 15;
theStimulusAreaDeg2 = theStimulusExtentDeg^2;

% Scale background to target cd/m2
%
% This makes units Watts/sr-m2-wlband
% Wavelength band is 2 here, which we need
% to keep track of.
targetLum = 1000;
theBGDeviceRawLum = T_xyz(2,:)*bgSpd;
theBgDeviceSpdScaled = targetLum*bgSpd/theBGDeviceRawLum;




