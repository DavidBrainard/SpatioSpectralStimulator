% SpectralTest
%
% Start exploring spectral fits with swubprimarys
%
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Define calibration filenames/params
projectorCalName = 'SACC';
projectorNInputLevels = 256;
primaryCalNames = {'SACCPrimary1' 'SACCPrimary1' 'SACCPrimary1'};
primaryNInputLevels = 252;

%% Load projector calibration
projectorCal = LoadCalFile(projectorCalName);
projectorCalObj = ObjectToHandleCalOrCalStruct(projectorCal);
CalibrateFitGamma(projectorCalObj, projectorNInputLevels);

%% Load subprimary calibrations
subprimaryCals = cell(3,1);
subprimaryCalObjs = cell(3,1);
for cc = 1:length(primaryCalNames)
    subprimaryCals{cc} = LoadCalFile(primaryCalNames{cc});
    subprimaryCalObjs{cc} = ObjectToHandleCalOrCalStruct(subprimaryCals{cc});
    CalibrateFitGamma(subprimaryCalObjs{cc}, primaryNInputLevels);
end

%% Get out some data to work with for now
S = subprimaryCalObjs{1}.get('S');
wls = SToWls(S);
ambient = subprimaryCalObjs{1}.get('P_ambient');
P_device = subprimaryCalObjs{1}.get('P_device');
gammaInput = subprimaryCalObjs{1}.get('gammaInput');
gammaTable = subprimaryCalObjs{1}.get('gammaTable');
gammaMeasurements = subprimaryCals{1}.rawData.gammaCurveMeanMeasurements;
[nPrimaries,nMeas,~] = size(gammaMeasurements);

%% Let's look at little at the subprimary calibration.
%
% Eventually this will be handled by the analyze program,
% when it is generalized for more than three primaries.  But
% we are impatient people so we will hack something up here.
PLOT_SUBPRIMARYINVARIANCE = false;
if (PLOT_SUBPRIMARYINVARIANCE)
    for pp = 1:nPrimaries
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

%% Plot subprimary gamma functions
PLOT_SUBPRIMARYGAMMA = false;
if (PLOT_SUBPRIMARYGAMMA)
    for pp = 1:nPrimaries
        figure; hold on;
        plot(subprimaryCals{1}.rawData.gammaInput,subprimaryCals{1}.rawData.gammaTable(:,pp),'ko','MarkerSize',12,'MarkerFaceColor','k');
        plot(gammaInput,gammaTable(:,pp),'k','LineWidth',2);
    end 
end

%% Background xy
targetBgxy = [0.3127 0.3290]';

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
ledContrastReMax = 0.025;
ledContrastReMaxWithHeadroom = 1.1*ledContrastReMax;
plotAxisLimit = 2;

% When we compute a specific image, we may not want full contrast available
% with the primaries. This tells us fraction of max available relative to
% ledContrastReMax.
imageModulationContrast = 0.0025/ledContrastReMax;

% Frozen factor for SRGB conversions, so it's preserved across contrasts
scaleFactor = 2.4985e+04;

% Image spatial parameters
sineFreq = 6;
gaborSdImageFraction = 0.1;

% Bit levels. We want to understand effects of bit depth.  Here we specify
% what bit depths we're using.
nSubprimaryLevels = 252;
displayBits = 8;
nDisplayLevels = 2^displayBits;

% This is a computational bit depth that we use to define the lookup table
% between contrast and primary values.
fineBits = 14;
nFineLevels = 2^fineBits;

% Image size
imageN = 512;

%% Generate a OneLight cal file with narrowband Subprimarys
% S = [380 2 201];
% cal = sssSpoofOneLightCal('S',S, ....
%     'plotBasis',false,...
%     'gaussianPeakWls',[437 485 540 562 585 618 652], ...
%     'gaussianFWHM',25);
% S = cal.describe.S;
% ambientSpd = cal.computed.pr650MeanDark;
% wls = SToWls(S);
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
nPrimaries = subprimaryCalObjs{1}.get('nDevices');
halfOnPrimaries = 0.5*ones(nPrimaries,1);
halfOnSpd = PrimaryToSpd(subprimaryCalObjs{1},halfOnPrimaries);

%% Make sure gamma correction behaves well
SetGammaMethod(subprimaryCalObjs{1},0);
halfOnSettings = PrimaryToSettings(subprimaryCalObjs{1},halfOnPrimaries);
halfOnPrimariesChk = SettingsToPrimary(subprimaryCalObjs{1},halfOnSettings);
if (max(abs(halfOnPrimaries-halfOnPrimariesChk)) > 1e-8)
    error('Gamma self-inversion not sufficiently precise');
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
bgSettings = round((nSubprimaryLevels-1)*bgPrimaries);
bgPrimaries = double(bgSettings)/(nSubprimaryLevels-1);
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
isolatingPrimaries1 = QuantizePrimaries(isolatingPrimaries1,(nSubprimaryLevels-1));

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
isolatingPrimaries2 = QuantizePrimaries(isolatingPrimaries2,(nSubprimaryLevels-1));

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
isolatingPrimaries3 = QuantizePrimaries(isolatingPrimaries3,(nSubprimaryLevels-1));

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
plot(wls,bgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd1,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd1,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd1(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd1(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 1');
ylim([0 2]);
subplot(2,2,3); hold on
plot(wls,bgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd2,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd2,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd2(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd2(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 2');
ylim([0 2]);
subplot(2,2,4); hold on
plot(wls,bgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd3,'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd3,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd3(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd3(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 3');
ylim([0 2]);

%% Create lookup table that maps [-1,1] to desired LMS contrast at a very fine scale
%
% Also find and save best mixture of quantized primaries to acheive those
% contrasts.
fprintf('Making fine contrast to LMS lookup table\n');
fineContrastLevels = linspace(-1,1,nFineLevels);
spdMatrix = [isolatingSpd1, isolatingSpd2, isolatingSpd3];
LMSMatrix = T_cones*spdMatrix;
for ll = 1:nFineLevels
    % Find the LMS values corresponding to desired contrast
    fineDesiredContrast(:,ll) = fineContrastLevels(ll)*ledContrastReMax*targetLMSContrast;
    fineDesiredLMS(:,ll) = ContrastToExcitation(fineDesiredContrast(:,ll),bgLMS);
    
    % Find primary mixture to best prodcue those values
    thisMixture = LMSMatrix\fineDesiredLMS(:,ll);
    thisMixture(thisMixture > 1) = 1;
    thisMixture(thisMixture < 0) = 0;
    
    % Store
    finePrimaries(:,ll) = thisMixture;
    finePredictedLMS(:,ll) = T_cones*spdMatrix*thisMixture;
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
finePtCloud = pointCloud(finePredictedLMS');
for ll = 1:nDisplayLevels
    quantizedLMSContrast(:,ll) = quantizedContrastLevels(ll)*ledContrastReMax*targetLMSContrast;
    quantizedLMS(:,ll) = ContrastToExcitation(quantizedLMSContrast(:,ll),bgLMS);
    
    minIndices(ll) = findNearestNeighbors(finePtCloud,quantizedLMS(:,ll)',1);
    predictedQuantizedLMS(:,ll) = finePredictedLMS(:,minIndices(ll));
    quantizedDisplayPrimaries(:,ll) = finePrimaries(:,minIndices(ll));
end

%% Make Gabor patch in range 0-1
%
% This is our contrast modulation
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawMonochromeSineImage = MakeSineImage(0,sineFreq,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawMonochromeGaborImage = imageModulationContrast*rawMonochromeSineImage.*gaussianWindow;

% Quantized for display bit depth
displayIntegerMonochromeGaborImage = PrimariesToIntegerPrimaries((rawMonochromeGaborImage+1)/2,nDisplayLevels);
displayIntegerMonochromeGaborCal = ImageToCalFormat(displayIntegerMonochromeGaborImage);

% Quantized for fine bit depth
fineIntegerMonochromeGaborImage = PrimariesToIntegerPrimaries((rawMonochromeGaborImage+1)/2,nFineLevels);
fineIntegerMonochromeGaborCal = ImageToCalFormat(fineIntegerMonochromeGaborImage);

%% Create the Gabor image with desired LMS contrasts
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


%% Create the Gabor image with quantized primary mixtures
fprintf('Making Gabor primary mixture image\n');
quantizedDisplayPrimariesGaborCal = zeros(3,imageN*imageN);
for ii = 1:imageN*imageN
    thisIndex = displayIntegerMonochromeGaborCal(ii);
    quantizedDisplayPrimariesGaborCal(:,ii) = quantizedDisplayPrimaries(:,thisIndex);
end

%% Convert of useful formats for analysis, rendering
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
%scaleFactor = max(quantizedSRGBPrimaryCal(:));
quantizedSRGBCal = SRGBGammaCorrect(quantizedSRGBPrimaryCal/scaleFactor,0);
quantizedSRGBImage = uint8(CalFormatToImage(quantizedSRGBCal,imageN,imageN));

% Show the SRGB image
figure; imshow(quantizedSRGBImage)

%% Plot slice through LMS contrast image
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

%% Light level tests
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




