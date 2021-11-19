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
% This is a standard calibration file for the DLP screen,
% with the subprimaries set to something.  As we'll see below,
% we're going to rewrite those.
screenCalName = 'SACC';
screenNInputLevels = 256;

% These are the calibration files for each of the primaries, which
% then entails measuring the spectra of all the subprimaries for that
% primary.
channelCalNames = {'SACCPrimary1' 'SACCPrimary1' 'SACCPrimary1'};
channelNInputLevels = 253;

%% Load screen calibration.
screenCal = LoadCalFile(screenCalName);
screenCalObj = ObjectToHandleCalOrCalStruct(screenCal);
CalibrateFitGamma(screenCalObj, screenNInputLevels);

%% Load channel calibrations.
nPrimaries = 3;
channelCals = cell(nPrimaries ,1);
channelCalObjs = cell(nPrimaries ,1);
for cc = 1:length(channelCalNames)
    channelCals{cc} = LoadCalFile(channelCalNames{cc});

    channelCalObjs{cc} = ObjectToHandleCalOrCalStruct(channelCals{cc});
    CalibrateFitGamma(channelCalObjs{cc}, channelNInputLevels);
end

%% Get out some data to work with.
%
% This is from the channel calibration file.
S = channelCalObjs{1}.get('S');
wls = SToWls(S);
ambientSpd = channelCalObjs{1}.get('P_ambient');
if (isempty(ambientSpd))
    channelCalObjs{1}.P_ambient = zeros(size(wls));
end
P_device = channelCalObjs{1}.get('P_device');
gammaInput = channelCalObjs{1}.get('gammaInput');
gammaTable = channelCalObjs{1}.get('gammaTable');
gammaMeasurements = channelCals{1}.rawData.gammaCurveMeanMeasurements;
[nChannels,nMeas,~] = size(gammaMeasurements);

%% Cone fundamentals and XYZ CMFs.
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
T_cones = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);
load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,S);

%% Let's look at little at the channel calibration.
%
% Eventually this will be handled by the analyze program,
% when it is generalized for more than three primaries.  But
% we are impatient people so we will hack something up here.
PLOT_CHANNELINVARIANCE = false;
if (PLOT_CHANNELINVARIANCE)
    for pp = 1:nChannels
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

%% Plot channel gamma functions.
PLOT_CHANNELGAMMA = false;
if (PLOT_CHANNELGAMMA)
    for pp = 1:nChannels
        figure; hold on;
        plot(channelCals{1}.rawData.gammaInput,channelCals{1}.rawData.gammaTable(:,pp),'ko','MarkerSize',12,'MarkerFaceColor','k');
        plot(gammaInput,gammaTable(:,pp),'k','LineWidth',2);
    end 
end

%% Plot x,y if desired.
PLOT_CHANNELCHROMATICITY = false;
if (PLOT_CHANNELCHROMATICITY)
    figure; hold on;
    for pp = 1:nChannels
        for mm = 1:nMeas
            % XYZ calculation for each measurement
            spd_temp = squeeze(gammaMeasurements(pp,mm,:));      
            XYZ_temp = T_xyz*spd_temp; 
            xyY_temp = XYZToxyY(XYZ_temp);
            
            plot(xyY_temp(1,:),xyY_temp(2,:),'r.','Markersize',10); % Coordinates of the channel
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
halfOnChannels = 0.5*ones(nChannels,1);
halfOnSpd = PrimaryToSpd(channelCalObjs{1},halfOnChannels);

%% Make sure gamma correction behaves well with unquantized conversion.
SetGammaMethod(channelCalObjs{1},0);
halfOnSettings = PrimaryToSettings(channelCalObjs{1},halfOnChannels);
halfOnPrimariesChk = SettingsToPrimary(channelCalObjs{1},halfOnSettings);
if (max(abs(halfOnChannels-halfOnPrimariesChk)) > 1e-8)
    error('Gamma self-inversion not sufficiently precise');
end

%% Use quantized conversion from here on.
%
% Comment in the line that refits the gamma to see
% effects of extreme quantization below
%
% CalibrateFitGamma(channelCalObjs{1},10);
SetGammaMethod(channelCalObjs{1},2);
SetGammaMethod(channelCalObjs{2},2);
SetGammaMethod(channelCalObjs{3},2);

%% Use extant machinery to get primaries from spectrum.
%
% This isn't used in our calculations.  Any difference in the
% two lines here reflects a bug in the SpdToPrimary/PrimaryToSpd pair.  
halfOnPrimariesChk = SpdToPrimary(channelCalObjs{1},halfOnSpd);
halfOnSpdChk = PrimaryToSpd(channelCalObjs{1},halfOnPrimariesChk);
figure; hold on;
plot(wls,halfOnSpd,'r','LineWidth',3);
plot(wls,halfOnSpdChk,'k','LineWidth',1);

%% Show effect of quantization.
%
% It's very small at the nominal 252 levels of the subprimaries, but will
% increase if you refit the gamma functios to a small number of levels.
halfOnPrimariesChk = SpdToPrimary(channelCalObjs{1},halfOnSpd);
halfOnSettingsChk = PrimaryToSettings(channelCalObjs{1},halfOnPrimariesChk);
halfOnPrimariesChk1 = SettingsToPrimary(channelCalObjs{1},halfOnSettingsChk);
halfOnSpdChk1 = PrimaryToSpd(channelCalObjs{1},halfOnPrimariesChk1);
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
    [bgPrimaries(:,pp),obtainedBgSpd(:,pp),obtainedBgXYZ(:,pp)] = FindBgChannelPrimaries(targetBgXYZ,T_xyz,channelCalObjs{pp}, ...
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
    [isolatingPrimaries(:,pp),isolatingPrimariesQuantized(:,pp),isolatingSpd(pp,:),isolatingContrast(pp,:)] = FindChannelPrimaries(targetMaxLMSContrast(:,pp), ...
        targetPrimaryHeadroom,targetContrastReMax,bgPrimaries(:,pp), ...
        T_cones,channelCalObjs{pp},B_natural,projectIndices,primaryHeadroom,targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
end


%% Measure the desired target primaries (THIS PART HAS BEEN UPDATED - SEMIN)
% This result will be used to compute the image.
%
% Make a loop for measuring all primaries.
for pp = 1:nPrimaries
    [targetScreenSpdMeasured(pp,:)] = MeasureChannelPrimaries(isolatingPrimaries(:,pp),channelNInputLevels,channelCalObjs{pp},pp,'projectorMode',true,'measurementOption',false,'verbose',false);
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

%% Set the screen primaries
%
% We want these to match those we set up with the
% channel calculations above.  Need to reset
% sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device',isolatingSpd');
SetSensorColorSpace(screenCalObj,T_cones,S)

%% Create lookup table that maps [-1,1] to desired LMS contrast at a very fine scale.
fprintf('Making fine contrast to LMS lookup table\n');
bgLMS = T_cones * sum(obtainedBgSpd,2);
fineContrastLevels = linspace(-1,1,nFineLevels);
fineDesiredContrast = targetContrastReMax*targetLMSContrast*fineContrastLevels;
fineDesiredLMS = ContrastToExcitation(fineDesiredContrast,bgLMS);
finePrimaries = PrimaryToGamut(screenCalObj,SensorToPrimary(screenCalObj,fineDesiredLMS));
finePredictedLMS = PrimaryToSensor(screenCalObj,finePrimaries);

% This is an older looped way of doing the above, which does not use
% the calbration infrastructure.  Delete when satisifed that the
% differences printout out are always small.
CHKNEWCODE = true;
if (CHKNEWCODE)
    spdMatrix = isolatingSpd';
    LMSMatrix = T_cones * isolatingSpd';
    for ll = 1:nFineLevels
        % Find primary mixture to best prodcue those values
        thisMixture = LMSMatrix\fineDesiredLMS(:,ll);
        thisMixture(thisMixture > 1) = 1;
        thisMixture(thisMixture < 0) = 0;

        % Store
        finePrimariesChk(:,ll) = thisMixture;
        finePredictedLMSChk(:,ll) = T_cones * spdMatrix * thisMixture;
    end
    if (max(abs(finePrimariesChk(:)-finePrimaries(:))) > 1e-12)
        error('Do not get same answer in two essentially the same ways');
    end
    if (max(abs(finePredictedLMSChk(:)-finePredictedLMS(:))) > 1e-13)
        error('Do not get same answer in two essentially the same ways');
    end  
end

%% DHB got to here in his quest to understand and update this code
% Do this at quantized levels
fprintf('Making display quantized primary lookup table\n');
quantizedIntegerLevels = 1:screenNInputLevels;
quantizedContrastLevels = (2*(quantizedIntegerLevels-1)/(screenNInputLevels-1))-1;
quantizedLMSContrast = zeros(3,screenNInputLevels);
quantizedLMS = zeros(3,screenNInputLevels);
minIndices = zeros(1,screenNInputLevels);
predictedQuantizedLMS = zeros(3,screenNInputLevels);
quantizedDisplayPrimaries = zeros(3,screenNInputLevels);

% Set up point cloud for fast finding of nearest neighbors
finePtCloud = pointCloud(finePredictedLMS');
for ll = 1:screenNInputLevels
    quantizedLMSContrast(:,ll) = quantizedContrastLevels(ll)*targetContrastReMax*targetLMSContrast;
    quantizedLMS(:,ll) = ContrastToExcitation(quantizedLMSContrast(:,ll),bgLMS);
    
    minIndices(ll) = findNearestNeighbors(finePtCloud,quantizedLMS(:,ll)',1);
    predictedQuantizedLMS(:,ll) = finePredictedLMS(:,minIndices(ll));
    quantizedDisplayPrimaries(:,ll) = finePrimaries(:,minIndices(ll));
end

%% Make Gabor patch in range 0-1.
%
% This is our contrast modulation
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawMonochromeSineImage = MakeSineImage(0,sineFreqCyclesPerImage,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawMonochromeGaborImage = imageModulationContrast*rawMonochromeSineImage.*gaussianWindow;

% Quantized for display bit depth
displayIntegerMonochromeGaborImage = PrimariesToIntegerPrimaries((rawMonochromeGaborImage+1)/2,screenNInputLevels);
displayIntegerMonochromeGaborCal = ImageToCalFormat(displayIntegerMonochromeGaborImage);

% Quantized for fine bit depth
fineIntegerMonochromeGaborImage = PrimariesToIntegerPrimaries((rawMonochromeGaborImage+1)/2,nFineLevels);
fineIntegerMonochromeGaborCal = ImageToCalFormat(fineIntegerMonochromeGaborImage);

%% Create the Gabor image with desired LMS contrasts.
fprintf('Making Gabor desired (fine) LMS contrast image\n');
quantizedFineLMSGaborCal = zeros(3,imageN*imageN);
for ii = 1:imageN*imageN
    thisIndex = fineIntegerMonochromeGaborImage(ii);
    fineLMSContrastCal(:,ii) = fineDesiredContrast(:,thisIndex);
    quantizedFineLMSGaborCal(:,ii) = finePredictedLMS(:,thisIndex);
end
fineLMSContrastGaborImage = CalFormatToImage(fineLMSContrastCal,imageN,imageN);
meanLMS = mean(quantizedFineLMSGaborCal,2);
quantizedFineContrastGaborCal = ExcitationsToContrast(quantizedFineLMSGaborCal,meanLMS);
quantizedFineContrastGaborImage = CalFormatToImage(quantizedFineContrastGaborCal,imageN,imageN);

%% Create the Gabor image with quantized primary mixtures.
fprintf('Making Gabor primary mixture image\n');
quantizedDisplayPrimariesGaborCal = zeros(3,imageN*imageN);
for ii = 1:imageN*imageN
    thisIndex = displayIntegerMonochromeGaborCal(ii);
    quantizedDisplayPrimariesGaborCal(:,ii) = quantizedDisplayPrimaries(:,thisIndex);
end

%% Convert of useful formats for analysis, rendering.
%
% Get spectral power distribution
fprintf('Convert Gabor for rendering, analysis\n');
quantizedSpdCal = spdMatrix*quantizedDisplayPrimariesGaborCal;

% Quantized LMS image and cone contrast image
quantizedLMSCal = T_cones*quantizedSpdCal;
meanLMS = mean(quantizedLMSCal,2);
quantizedContrastCal = ExcitationsToContrast(quantizedLMSCal,meanLMS);
quantizedContrastImage = CalFormatToImage(quantizedContrastCal,imageN,imageN);

% SRGB image via XYZ
quantizedXYZCal = T_xyz*quantizedSpdCal;
quantizedSRGBPrimaryCal = XYZToSRGBPrimary(quantizedXYZCal);
scaleFactor = max(quantizedSRGBPrimaryCal(:));
quantizedSRGBCal = SRGBGammaCorrect(quantizedSRGBPrimaryCal/(2*scaleFactor),0);
quantizedSRGBImage = uint8(CalFormatToImage(quantizedSRGBCal,imageN,imageN));

% Show the SRGB image
figure; imshow(quantizedSRGBImage)

%% Now compute screen image.
%
% First step is to make a DLP calibration file that has as primaries
% the three spds we've computed above.
%
% In an actual display program, we would set each of the primary's
% subprimaries to isolatingPrimaries1, isolatingPrimaries2,
% isolatingPrimaries3 as computed above.  That now allows the DLP
% to produce mixtures of these primaries.  Here we tell the calibration
% object for the DLP that it has these desired primaries.
P_device = isolatingSpd';
screenCal.processedData.P_device = P_device;

% Initialze the calibration structure
screenCal = SetGammaMethod(screenCal,2);

% Convert excitations image to screen settings
[screenSettingsCal,outOfGamutIndex] = SensorToSettings(screenCal,quantizedFineLMSGaborCal);
if (any(outOfGamutIndex))
    error('Oops: Some pixels out of gamut');
end
screenSettingsImage = CalFormatToImage(screenSettingsCal,imageN,imageN);
figure; clf;
imshow(screenSettingsImage)

% Show this image on the DLP, and it should look more or less like
% the sRGB image we display below.
testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
testFilename = fullfile(testFiledir,'testImageData1');
save(testFilename,'screenSettingsImage','isolatingPrimaries');

%% Plot slice through LMS contrast image.
figure; hold on
plot(1:imageN,100*quantizedContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*fineLMSContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);
%plot(1:imageN,100*quantizedFineContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*quantizedContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*fineLMSContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);
%plot(1:imageN,100*quantizedFineContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*quantizedContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*fineLMSContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
%plot(1:imageN,100*quantizedFineContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Measure LMS contrast on the gabor patch image. (THIS PART HAS BEEN ADDED - SEMIN)
[targetLMSContrastMeasured] = MeasureLMSContrastGaborPatch(quantizedContrastImage,isolatingPrimaries,screenCalObj,obtainedBgSpd,T_cones, ...
    channelNInputLevels,channelCalObjs,'projectorMode',true,'measurementOption',false,'verbose',true);

%% DAVID - Add plot of primaries.

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




