% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Parameters

%% Background xy
targetBgxy = [0.325 0.34]';

% Target color direction and max contrasts.
%
% This is the basic desired modulation direction positive excursion. We go
% equally in positive and negative directions.
targetLMSContrast = [1 -1 0]';

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
ledContrastReMax = 0.05;
ledContrastReMaxWithHeadroom = 1.1*ledContrastReMax;

% When we compute a specific image, we may not want full contrast available
% with the primaries. This tells us fraction of max available relative to
% ledContrastReMax.
imageModulationContrast = 0.25;

% Image spatial parameters
sineFreq = 6;
gaborSdImageFraction = 0.15;

% Bit levels. We want to understand effects of bit depth.  Here we specify
% what bit depths we're using.
LEDBits = 8;
nLEDLevels = 2^LEDBits;
displayBits = 8;
nDisplayLevels = 2^displayBits;

% This is a computational bit depth that we use to define the lookup table
% between contrast and primary values.
fineBits = 14;
nFineLevels = 2^fineBits;

% Image size
imageN = 512;

%% Generate a OneLight cal file with narrowband LEDs
S = [380 2 201];
cal = sssSpoofOneLightCal('S',S, ....
    'plotBasis',false,...
    'gaussianPeakWls',[437 485 540 562 585 618 652], ...
    'gaussianFWHM',25);
S = cal.describe.S;
ambientSpd = cal.computed.pr650MeanDark;
wls = SToWls(S);
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(wls > lowProjectWl & wls < highProjectWl);

%% Cone fundamentals and XYZ CMFs
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
T_cones = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);
load T_xyz1931
T_xyz = 683*SplineCmf(S_xyz1931,T_xyz1931,S);

%% Get half on spectrum
%
% This is useful for scaling things reasonably
halfOnPrimaries = 0.5*ones(cal.describe.numWavelengthBands,1);
halfOnSpd = OLPrimaryToSpd(cal,halfOnPrimaries);
nPrimaries = size(halfOnPrimaries,1);

%% Make sure gamma correction behaves well
halfOnSettings = OLPrimaryToSettings(cal,halfOnPrimaries);
if (max(abs(halfOnPrimaries-halfOnSettings)) > 1e-8)
    error('Spoofed cal file does not have linear gamma');
end
halfOnPrimariesChk = OLSettingsToPrimary(cal,halfOnSettings);
if (max(abs(halfOnPrimaries-halfOnPrimariesChk)) > 1e-8)
    error('Spoofed cal file does not have linear inverse gamma');
end

%% Use OL machinery to get primaries from spectrum
halfOnPrimariesChk = OLSpdToPrimary(cal,halfOnSpd);
halfOnSpdChk = OLPrimaryToSpd(cal,halfOnPrimariesChk);
% figure; hold on;
% plot(wls,halfOnSpd,'r','LineWidth',3);
% plot(wls,halfOnSpdChk,'k','LineWidth',1);

%% Set up basis to try to keep spectra close to
basisType = 'fourier';
switch (basisType)
    case 'cieday'
        load B_cieday
        B_natural = SplineSpd(S_cieday,B_cieday,S);
    case 'fourier'
        B_natural = MakeFourierBasis(S,7);
    otherwise
        error('Unknown basis set specified');
end

%% Some transformation matrices
M_NaturalToXYZ = T_xyz*B_natural(:,1:3);
M_XYZToNatural = inv(M_NaturalToXYZ);
M_NaturalToLMS = T_cones*B_natural(:,1:3);
M_LMSToNatural = inv(M_NaturalToLMS);

%% Get background
targetBgxyY = [targetBgxy ; 1];
targetBgXYZ = xyYToXYZ(targetBgxyY);
desiredBgSpd = B_natural(:,1:3)*M_XYZToNatural*targetBgXYZ;
desiredBgSpd = sum(halfOnSpd)*desiredBgSpd/sum(desiredBgSpd);
desiredBgXYZ = T_xyz*desiredBgSpd;
desiredBgxyY = XYZToxyY(desiredBgXYZ);

%% Search for desired background
B_primary = cal.computed.pr650M;
startingBgPrimaries = OLSpdToPrimary(cal,desiredBgSpd);
startingBgSpd = OLPrimaryToSpd(cal,startingBgPrimaries);
startingBgXYZ = T_xyz*startingBgSpd;
startingBgxyY = XYZToxyY(startingBgXYZ);

primaryHeadRoom = 0;
targetLambda = 10;
[bgPrimariesIncr] = ReceptorIsolateSpectral(T_xyz,desiredBgXYZ,B_primary,startingBgPrimaries,startingBgPrimaries, ...
    primaryHeadRoom,B_natural,projectIndices,targetLambda,ambientSpd,'EXCITATIONS',true);
bgPrimaries = startingBgPrimaries + bgPrimariesIncr;
bgSettings = round((nLEDLevels-1)*bgPrimaries);
bgPrimaries = double(bgSettings)/(nLEDLevels-1);
bgSpd = OLPrimaryToSpd(cal,bgPrimaries);
bgLMS = T_cones*bgSpd;
bgXYZ = T_xyz*bgSpd;
bgxyY = XYZToxyY(bgXYZ);
fprintf('Desired  background x,y = %0.3f,%0.3f\n',desiredBgxyY(1),desiredBgxyY(2));
fprintf('Starting background x,y = %0.3f,%0.3f\n',startingBgxyY(1),startingBgxyY(2));
fprintf('Obtained background x,y = %0.3f,%0.3f\n',bgxyY(1),bgxyY(2));
fprintf('Mean value of background primaries: %0.2f\n',mean(bgPrimaries));

%% Get primaries based on contrast specification
%
LMSContrast1 = ledContrastReMaxWithHeadroom*target1MaxLMSContrast;
targetLambda = 10;
[isolatingModulationPrimaries1] = ReceptorIsolateSpectral(T_cones,LMSContrast1,B_primary,bgPrimaries,bgPrimaries, ...
    primaryHeadRoom,B_natural,projectIndices,targetLambda,ambientSpd);
isolatingPrimaries1 = isolatingModulationPrimaries1 + bgPrimaries;

% Quantize
isolatingSettings1 = round((nLEDLevels-1)*isolatingPrimaries1);
isolatingPrimaries1 = double(isolatingSettings1)/(nLEDLevels-1);

% Report
isolatingSpd1 = OLPrimaryToSpd(cal,isolatingPrimaries1);
isolatingLMS1 = T_cones*isolatingSpd1;
isolatingContrast1 = ExcitationsToContrast(isolatingLMS1,bgLMS);
fprintf('Desired/obtained contrasts 1\n');
for rr = 1:length(LMSContrast1)
    fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,LMSContrast1(rr),isolatingContrast1(rr));
end
fprintf('Min/max primaries 1: %0.4f, %0.4f\n', ...
    min(isolatingPrimaries1), max(isolatingPrimaries1));

% Primary 2
LMSContrast2 = ledContrastReMaxWithHeadroom*target2MaxLMSContrast;
targetLambda = 10;
[isolatingModulationPrimaries2] = ReceptorIsolateSpectral(T_cones,LMSContrast2,B_primary,bgPrimaries,bgPrimaries, ...
    primaryHeadRoom,B_natural,projectIndices,targetLambda,ambientSpd);
isolatingPrimaries2 = isolatingModulationPrimaries2 + bgPrimaries;

% Quantize
isolatingSettings2 = round((nLEDLevels-1)*isolatingPrimaries2);
isolatingPrimaries2 = double(isolatingSettings2)/(nLEDLevels-1);

% Report
isolatingSpd2 = OLPrimaryToSpd(cal,isolatingPrimaries2);
isolatingLMS2 = T_cones*isolatingSpd2;
isolatingContrast2 = ExcitationsToContrast(isolatingLMS2,bgLMS);
fprintf('Desired/obtained contrasts 2\n');
for rr = 1:length(LMSContrast2)
    fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,LMSContrast2(rr),isolatingContrast2(rr));
end
fprintf('Min/max primaries 2: %0.4f, %0.4f\n', ...
    min(isolatingPrimaries2), max(isolatingPrimaries2));

% Primary 3
LMSContrast3 = ledContrastReMaxWithHeadroom*target3MaxLMSContrast;
targetLambda = 10;
[isolatingModulationPrimaries3] = ReceptorIsolateSpectral(T_cones,LMSContrast3,B_primary,bgPrimaries,bgPrimaries, ...
    primaryHeadRoom,B_natural,projectIndices,targetLambda,ambientSpd);
isolatingPrimaries3 = isolatingModulationPrimaries3 + bgPrimaries;

% Quantize
isolatingSettings3 = round((nLEDLevels-1)*isolatingPrimaries3);
isolatingPrimaries3 = double(isolatingSettings3)/(nLEDLevels-1);

% Report
isolatingSpd3 = OLPrimaryToSpd(cal,isolatingPrimaries3);
isolatingLMS3 = T_cones*isolatingSpd3;
isolatingContrast3 = ExcitationsToContrast(isolatingLMS3,bgLMS);
fprintf('Desired/obtained contrasts 3\n');
for rr = 1:length(LMSContrast3)
    fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,LMSContrast3(rr),isolatingContrast3(rr));
end
fprintf('Min/max primaries 3: %0.4f, %0.4f\n', ...
    min(isolatingPrimaries3), max(isolatingPrimaries3));


%% How close are spectra to subspace defined by basis?
theBgNaturalApproxSpd = B_natural*(B_natural(projectIndices,:)\bgSpd(projectIndices));
isolatingNaturalApproxSpd1 = B_natural*(B_natural(projectIndices,:)\isolatingSpd1(projectIndices));
isolatingNaturalApproxSpd2 = B_natural*(B_natural(projectIndices,:)\isolatingSpd2(projectIndices));
isolatingNaturalApproxSpd3 = B_natural*(B_natural(projectIndices,:)\isolatingSpd3(projectIndices));

% Plot
figure; clf; 
subplot(2,2,1); hold on
plot(wls,bgSpd,'b','LineWidth',2);
plot(wls,theBgNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),bgSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theBgNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Background');
ylim([0 2]);
subplot(2,2,2); hold on
plot(wls,isolatingSpd1,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd1,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd1(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd1(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 1');
ylim([0 2]);
subplot(2,2,3); hold on
plot(wls,isolatingSpd2,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd2,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd2(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd2(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 2');
ylim([0 2]);
subplot(2,2,4); hold on
plot(wls,isolatingSpd3,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd3,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd3(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd3(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 3');
ylim([0 2]);

%% Create lookup table that maps between [-1,1] in contrast to LMS
fprintf('Making fine contrast to LMS lookup table\n');
fineContrastLevels = linspace(-1,1,nFineLevels);
% spdMatrix1 = [bgSpd , theIsolatingDeviceSpdLower];
% spdMatrix2 = [bgSpd , isolatingSpd1];
spdMatrix = [isolatingSpd1, isolatingSpd2, isolatingSpd3];
LMSMatrix = T_cones*spdMatrix;
for ll = 1:nFineLevels
    % Find the LMS values corresponding to desired contrast
    fineLMSContrast(:,ll) = fineContrastLevels(ll)*ledContrastReMax*targetLMSContrast;
    fineLMS(:,ll) = ContrastToExcitation(fineLMSContrast(:,ll),bgLMS);
    
    % Find primary mixture to best prodcue those values
    thisMixture = LMSMatrix\fineLMS(:,ll);
    thisMixture(thisMixture > 1) = 1;
    thisMixture(thisMixture < 0) = 0;
    
    % Store
    finePrimaries(:,ll) = thisMixture;
    predictedFineLMS(:,ll) = T_cones*spdMatrix*thisMixture;
end

% Do this at quantized levels
fprintf('Making display quantized primary lookup table\n');
quantizedIntegerLevels = 1:nDisplayLevels;
quantizedContrastLevels = (2*(quantizedIntegerLevels-1)/(nDisplayLevels-1))-1;
quantizedLMSContrast = zeros(3,nDisplayLevels);
quantizedLMS = zeros(3,nDisplayLevels);
minIndices = zeros(1,nDisplayLevels);
predictedQuantizedLMS = zeros(3,nDisplayLevels);
quantizedDisplayPrimaries = zeros(3,nDisplayLevels);

% Set up point cloud for fast finding of nearest neighbors
finePtCloud = pointCloud(fineLMS');
for ll = 1:nDisplayLevels
    quantizedLMSContrast(:,ll) = quantizedContrastLevels(ll)*ledContrastReMax*targetLMSContrast;
    quantizedLMS(:,ll) = ContrastToExcitation(quantizedLMSContrast(:,ll),bgLMS);
    
    minIndices(ll) = findNearestNeighbors(finePtCloud,quantizedLMS(:,ll)',1);
    predictedQuantizedLMS(:,ll) = fineLMS(:,minIndices(ll));
    quantizedDisplayPrimaries(:,ll) = finePrimaries(:,minIndices(ll));      
end

%% Make Gabor patch in range 0-1
%
% This is our contrast modulation
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawSineImage = MakeSineImage(0,sineFreq,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawGaborImage = imageModulationContrast*rawSineImage.*gaussianWindow;
quantizedIntegerGaborImage = round((nDisplayLevels-1)*(rawGaborImage+1)/2 );
quantizedIntegerGaborImageCal = ImageToCalFormat(quantizedIntegerGaborImage);

%% Create the Gabor image with quantized primary mixtures
fprintf('Making Gabor primary mixture image\n');
quantizedDisplayPrimariesGaborImageCal = zeros(3,imageN*imageN);
for ii = 1:imageN*imageN
    thisIndex = quantizedIntegerGaborImageCal(ii);
    quantizedDisplayPrimariesGaborImageCal(:,ii) = quantizedDisplayPrimaries(:,thisIndex);
end

%% Convert of useful formats for analysis, rendering
%
% Get spectral power distribution
fprintf('Convert Gabor for rendering, analysis\n');
isolatingSpdCal = spdMatrix*quantizedDisplayPrimariesGaborImageCal;

% LMS image and cone contrast image
isolatingLMSCal = T_cones*isolatingSpdCal;
isolatingLMSImage = CalFormatToImage(isolatingLMSCal,imageN,imageN);
meanLMS = mean(isolatingLMSCal,2);
isolatingContrastCal = zeros(3,imageN*imageN);
for ii = 1:imageN*imageN
    isolatingContrastCal(:,ii) = ExcitationsToContrast(isolatingLMSCal(:,ii),meanLMS);
end
isolatingContrastImage = CalFormatToImage(isolatingContrastCal,imageN,imageN);

% SRGB image via XYZ
isolatingXYZCal = T_xyz*isolatingSpdCal;
isolatingSRGBPrimaryCal = XYZToSRGBPrimary(isolatingXYZCal);
scaleFactor = max(isolatingSRGBPrimaryCal(:));
isolatingSRGBCal = SRGBGammaCorrect(isolatingSRGBPrimaryCal/scaleFactor,0);
isolatingSRGBImage = uint8(CalFormatToImage(isolatingSRGBCal,imageN,imageN));

% Show the SRGB image
figure; imshow(isolatingSRGBImage)

%% Plot slice through LMS contrast image
figure; hold on
plot(1:imageN,100*isolatingContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*isolatingContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*isolatingContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
title('Image Slice, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');

%% Light level tests

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




