% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Parameters
targetMaxLMSContrast = [1 -1 0]'/10;
primaryContrastReMax = 0.1;
imageModulationContrast = 0.5;
sineFreq = 4;
gaborSdImageFraction = 0.1;

% Bit levels
LEDBits = 10;
nLEDLevels = 2^LEDBits;
displayBits = 8;
nDisplayLevels = 2^displayBits;
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
halfOnPrimaries = 0.5*ones(cal.describe.numWavelengthBands,1);
halfOnSpd = OLPrimaryToSpd(cal,halfOnPrimaries);

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

%% Fit a spectrum with the basis
%
% Generate a black body radiator and put it into right general scale
theTemp = 4500;
theBackgroundNaturalSpd1 = GenerateBlackBody(theTemp,wls);
theTemp2 = 6000;
theBackgroundNaturalSpd2 = GenerateBlackBody(theTemp2,wls);
theBackgroundNaturalSpd = theBackgroundNaturalSpd1 + theBackgroundNaturalSpd2;
theBackgroundNaturalSpd = ones(size(theBackgroundNaturalSpd));
theBackgroundNaturalSpd = sum(halfOnSpd)*theBackgroundNaturalSpd/sum(theBackgroundNaturalSpd);

% Use OL machinery to get primaries and then the spd
theIsolatingNaturalPrimaries = OLSpdToPrimary(cal,theBackgroundNaturalSpd,'lambda',0.01);
theFitSpd = OLPrimaryToSpd(cal,theIsolatingNaturalPrimaries);

% Use non-neg least squares for comparison
weights = lsqnonneg(cal.computed.pr650M,theBackgroundNaturalSpd);
theFitSpdRegress = cal.computed.pr650M*weights;

% Plot
% figure; clf; hold on
% plot(wls,theIsolatingNaturalSpd,'k','LineWidth',3);
% plot(wls,theFitSpd,'r','LineWidth',2);
% plot(wls,theFitSpdRegress,'g','LineWidth',1);

%% Generate cone isolating spectral modulations with a simple basis
basisType = 'fourier';
switch (basisType)
    case 'cieday'
        load B_cieday
        theNaturalBasis = SplineSpd(S_cieday,B_cieday,S);
    case 'fourier'
        theNaturalBasis = MakeFourierBasis(S,7);
    otherwise
        error('Unknown basis set specified');
end

% Transform
M_NaturalPrimariesToCones = T_cones*theNaturalBasis(:,1:size(T_cones,1));
M_ConesToNaturalPrimaries = inv(M_NaturalPrimariesToCones);

% Get background within basis
theBgNaturalPrimaries = theNaturalBasis\theBackgroundNaturalSpd;
theBgNaturalSpd = theNaturalBasis*theBgNaturalPrimaries;
theBgNaturalLmsTarget = T_cones*theBackgroundNaturalSpd;
theBgNaturalXYZTarget = T_xyz*theBackgroundNaturalSpd;
figure; clf; hold on
plot(wls,theBgNaturalSpd,'b:','LineWidth',2);

%% Approximation to desired background
IDENTITY = false;
if (IDENTITY)
    B_DevicePrimary = eye(size(cal.computed.pr650M,1));
    theBgDevicePrimariesApprox = theBgNaturalSpd;
    initialDevicePrimariesApprox = theBackgroundNaturalSpd;
else
    B_DevicePrimary = cal.computed.pr650M;
    theBgDevicePrimariesApprox = OLSpdToPrimary(cal,theBgNaturalSpd);
    initialDevicePrimariesApprox = theBgDevicePrimariesApprox;
end
fprintf('Mean value of background primaries: %0.2f\n',mean(theBgDevicePrimariesApprox));

%% Search for desired background
ambientSpd = cal.computed.pr650MeanDark;
primaryHeadRoom = 0;
targetBasis = theNaturalBasis;
targetLambda = 10;
[backgroundDevicePrimariesIncr] = ReceptorIsolateSpectral(T_cones,theBgNaturalLmsTarget',B_DevicePrimary,theBgDevicePrimariesApprox,initialDevicePrimariesApprox, ...
    primaryHeadRoom,targetBasis,projectIndices,targetLambda,ambientSpd,'EXCITATIONS',true);
theBgDevicePrimaries = theBgDevicePrimariesApprox + backgroundDevicePrimariesIncr;
theBgDeviceSettings = round((nLEDLevels-1)*theBgDevicePrimaries);
theBgDevicePrimaries = double(theBgDeviceSettings)/(nLEDLevels-1);
theBgDeviceSpd = OLPrimaryToSpd(cal,theBgDevicePrimaries);
theBgDeviceLMS = T_cones*theBgDeviceSpd;
theBgDeviceXYZ = T_xyz*theBgDeviceSpd;

fprintf('Desired/obtained background excitations\n');
for rr = 1:length(theBgNaturalLmsTarget)
    fprintf('\tReceptor %d (desired/obtained): %0.2f, %0.2f\n',rr,theBgNaturalLmsTarget(rr),theBgDeviceLMS(rr));
end

% Define desired cone contrast
theLMSContrast = primaryContrastReMax*targetMaxLMSContrast;
targetLambda = 10;
initialDevicePrimaries = theBgDevicePrimaries;
[isolatingModulationDevicePrimaries] = ReceptorIsolateSpectral(T_cones,theLMSContrast',B_DevicePrimary,theBgDevicePrimaries,initialDevicePrimaries, ...
    primaryHeadRoom,targetBasis,projectIndices,targetLambda,ambientSpd);

isolatingDevicePrimariesUpper = isolatingModulationDevicePrimaries + theBgDevicePrimaries;
isolatingDeviceSettingsUpper = round((nLEDLevels-1)*isolatingDevicePrimariesUpper);
isolatingDevicePrimariesUpper = double(isolatingDeviceSettingsUpper)/(nLEDLevels-1);

isolatingDevicePrimariesLower = -isolatingModulationDevicePrimaries + theBgDevicePrimaries;
isolatingDeviceSettingsLower = round((nLEDLevels-1)*isolatingDevicePrimariesLower);
isolatingDevicePrimariesLower = double(isolatingDeviceSettingsLower)/(nLEDLevels-1);

theIsolatingDeviceSpdUpper = OLPrimaryToSpd(cal,isolatingDevicePrimariesUpper);
theIsolatingDeviceSpdLower = OLPrimaryToSpd(cal,isolatingDevicePrimariesLower);

theIsolatingDeviceLMS = T_cones*theIsolatingDeviceSpdUpper;
theIsolatingContrast = ExcitationsToContrast(theIsolatingDeviceLMS,theBgDeviceLMS);
fprintf('Desired/obtained contrasts\n');
for rr = 1:length(theLMSContrast)
    fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,theLMSContrast(rr),theIsolatingContrast(rr));
end
fprintf('Min/max primaries upper: %0.4f, %0.4f, lower: %0.4f, %0.4f\n', ...
    min(isolatingDevicePrimariesUpper), max(isolatingDevicePrimariesUpper), ...
    min(isolatingDevicePrimariesLower),  max(isolatingDevicePrimariesLower));

% How close are spectra to subspace defined by basis?
theBgNaturalApproxSpd = targetBasis*(targetBasis(projectIndices,:)\theBgDeviceSpd(projectIndices));
theIsolatingNaturalApproxSpdUpper = targetBasis*(targetBasis(projectIndices,:)\theIsolatingDeviceSpdUpper(projectIndices));
theIsolatingNaturalApproxSpdLower = targetBasis*(targetBasis(projectIndices,:)\theIsolatingDeviceSpdLower(projectIndices));

figure; clf; 
subplot(1,3,1); hold on
plot(wls,theBgDeviceSpd,'b','LineWidth',2);
plot(wls,theBgNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),theBgDeviceSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theBgNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Background');
ylim([0 2]);
subplot(1,3,2); hold on
plot(wls,theIsolatingDeviceSpdUpper,'b','LineWidth',2);
plot(wls,theIsolatingNaturalApproxSpdUpper,'r:','LineWidth',1);
plot(wls(projectIndices),theIsolatingDeviceSpdUpper(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theIsolatingNaturalApproxSpdUpper(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Modulated (+)');
ylim([0 2]);
subplot(1,3,3); hold on
plot(wls,theIsolatingDeviceSpdLower,'b','LineWidth',2);
plot(wls,theIsolatingNaturalApproxSpdLower,'r:','LineWidth',1);
plot(wls(projectIndices),theIsolatingDeviceSpdLower(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theIsolatingNaturalApproxSpdLower(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Modulated (-)');
ylim([0 2]);

%% Render quantized colors between the two primaries
centerN = imageN/2;

% Get the quantized primaries
nPrimaries = size(isolatingDevicePrimariesUpper,1);
for ii = 1:nPrimaries
    thesePrimaries = linspace(isolatingDevicePrimariesLower(ii),isolatingDevicePrimariesUpper(ii),256);
    for jj = 1:length(thesePrimaries)
        isolatingPrimariesRaw(ii,jj) = thesePrimaries(jj);
    end
end
isolatingSettings = uint8(255*isolatingPrimariesRaw);
isolatingPrimaries = double(isolatingSettings)/255;

% Get the quantized spds and LMS values
for ii = 1:size(isolatingPrimaries,2)
    isolatingSpd(:,ii) = OLPrimaryToSpd(cal,isolatingPrimaries(:,ii));
end
isolatingLMS = T_cones*isolatingSpd;

% Get the quantized XYZ values and convert to sRGB
isolatingXYZ = T_xyz*isolatingSpd;
for ii = 1:size(isolatingPrimaries,2)
    isolatingSRGBPrimary(:,ii) = XYZToSRGBPrimary(isolatingXYZ(:,ii));
end
scaleFactor = max(isolatingSRGBPrimary(:));
for ii = 1:size(isolatingPrimaries,2)
    isolatingSRGB(:,ii) = SRGBGammaCorrect(isolatingSRGBPrimary(:,ii)/scaleFactor,0);
end

%% Make Gabor patch
gaborSdPixels = gaborSdImageFraction*imageN;
rawSineImage = MakeSineImage(0,sineFreq,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawGaborImage = imageModulationContrast*rawSineImage.*gaussianWindow;
quantizedIntegerGaborImage = round((nDisplayLevels-1)*(rawGaborImage+1)/2 );
quantizedContrastGaborImage = (2*(quantizedIntegerGaborImage-1)/(nDisplayLevels-1))-1;

%% Create lookup table that maps between [-1,1] in contrast to LMS
fineContrastLevels = linspace(-1,1,nFineLevels);
spdMatrix1 = [theBgDeviceSpd , theIsolatingDeviceSpdLower];
spdMatrix2 = [theBgDeviceSpd , theIsolatingDeviceSpdUpper];
spdMatrix = [theBgDeviceSpd , theIsolatingDeviceSpdLower , theIsolatingDeviceSpdUpper];
LMSMatrix1 = T_cones*spdMatrix1;
LMSMatrix2 = T_cones*spdMatrix2;
LMSMatrix = T_cones*spdMatrix;
for ll = 1:nFineLevels
    fineLMSContrast(:,ll) = fineContrastLevels(ll)*theLMSContrast;
    fineLMS(:,ll) = ContrastToExcitation(fineLMSContrast(:,ll),theBgDeviceLMS);
    
    thisMixtureRaw1 = LMSMatrix1\fineLMS(:,ll);
    thisMixtureRaw2 = LMSMatrix2\fineLMS(:,ll);
    if (any(thisMixtureRaw1 < 0))
        thisMixture = [thisMixtureRaw2(1) 0 thisMixtureRaw2(2)]';
    else
        thisMixture = [thisMixtureRaw1(1) thisMixtureRaw1(2) 0]';
    end
    thisMixture(thisMixture > 1) = 1;
    thisMixture(thisMixture < 0) = 0;
    finePrimaries(:,ll) = thisMixture;
    predictedFineLMS(:,ll) = T_cones*spdMatrix*thisMixture;
end

% Do this at quantized levels
quantizedIntegerLevels = 1:nDisplayLevels;
quantizedContrastLevels = (2*(quantizedIntegerLevels-1)/(nDisplayLevels-1))-1;
for ll = 1:nDisplayLevels
    quantizedLMSContrast(:,ll) = quantizedContrastLevels(ll)*theLMSContrast;
    quantizedLMS(:,ll) = ContrastToExcitation(quantizedLMSContrast(:,ll),theBgDeviceLMS);
    minErr = Inf;
    minIndex = Inf;
    for jj = 1:nFineLevels
        thisDiff = quantizedLMS(:,ll)-fineLMS(:,jj);
        thisErr = norm(thisDiff);
        if (thisErr < minErr)
            minErr = thisErr;
            minIndex = jj;
        end
    end
    minIndices(ll) = minIndex;
    predictedQuantizedLMS(:,ll) = fineLMS(:,minIndex);
    quantizedPrimaries(:,ll) = finePrimaries(:,minIndex);      
end


% %% Create the continuous LMS image that we want
% for ii = 1:imageN
%     for jj = 1:imageN
%         for cc = 1:3
%             gaborImageLMS(ii,jj,cc) = theBgDeviceLMS(cc)*theLMSContrast(cc)*rawGaborImage(ii,jj) + theBgDeviceLMS(cc);
%         end
%     end
% end

%% Create the Gabor image with quantized primary mixtures
quantizedPrimariesGaborImage = zeros(imageN,imageN,3);
for ii = 1:imageN
    for jj = 1:imageN
        thisIndex = quantizedIntegerGaborImage(ii,jj);
        thisPrimaries = quantizedPrimaries(:,thisIndex);
        quantizedPrimariesGaborImage(ii,jj,:) = thisPrimaries;
    end
end

% gaborImagePrimariesRaw = zeros(imageN,imageN,3);
% fineLMS = zeros(3,imageN^2);
% predictLMS = zeros(3,imageN^2);
% inIndex = 1;
% for ii = 1:imageN
%     for jj = 1:imageN
%         fineLMS(:,inIndex) = squeeze(gaborImageLMS(ii,jj,:));
%         thisMixtureRaw1 = LMSMatrix1\fineLMS(:,inIndex);
%         thisMixtureRaw2 = LMSMatrix2\fineLMS(:,inIndex);
%         if (any(thisMixtureRaw1 < 0))
%             thisMixture = [thisMixtureRaw2(1) 0 thisMixtureRaw2(2)]';
%         else
%             thisMixture = [thisMixtureRaw1(1) thisMixtureRaw1(2) 0]';
%         end
%         thisMixture(thisMixture > 1) = 1;
%         thisMixture(thisMixture < 0) = 0;
%         predictLMS(:,inIndex) = T_cones*spdMatrix*thisMixture;
%         finePrimaries(:,inIndex) = thisMixture;
%         gaborImagePrimariesRaw(ii,jj,:) = thisMixture;
%         inIndex = inIndex+1;
%     end
% end
% gaborImagePrimaries = round((nDisplayLevels-1)*gaborImagePrimariesRaw);
% gaborImagePrimaries01 = gaborImagePrimaries/(nDisplayLevels-1);

% Convert gabor to SRGB, XYZ and LMS images
isolatingSRGBPrimaryImage = zeros(imageN,imageN,3);
isolatingXYZImage = zeros(imageN,imageN,3);
isolatingLMSImage = zeros(imageN,imageN,3);
predictLMSQuantized = zeros(3,imageN^2);
inIndex = 1;
for ii = 1:imageN
    for jj = 1:imageN
        thesePrimaries = squeeze(quantizedPrimariesGaborImage(ii,jj,:));
        quantizedPrimaries(:,inIndex) = thesePrimaries;
        thisSpd = spdMatrix*thesePrimaries;
        thisLMS = T_cones*thisSpd;
        thisXYZ = T_xyz*thisSpd;
        thisSRGBPrimary = XYZToSRGBPrimary(thisXYZ);
        isolatingXYZImage(ii,jj,:) = thisXYZ;
        isolatingLMSImage(ii,jj,:) = thisLMS;
        isolatingSRGBPrimaryImage(ii,jj,:) = thisSRGBPrimary;
        predictLMSQuantized(:,inIndex) = thisLMS;
        inIndex = inIndex + 1;
    end
end
scaleFactor = max(isolatingSRGBPrimaryImage(:));
isolatingSRGB = uint8(zeros(imageN,imageN,3));
for ii = 1:imageN
    for jj = 1:imageN
        isolatingSRGB(ii,jj,:) = SRGBGammaCorrect(isolatingSRGBPrimaryImage(ii,jj,:)/scaleFactor,0);
    end
end

%% Compute some errors
% Lerror = predictLMS(1,:)-predictLMSQuantized(1,:);
% quantError = finePrimaries(1,:) - quantizedPrimaries(1,:);
% figure; plot((nDisplayLevels-1)*quantError,Lerror,'ro');


% Show the SRGB image
figure; imshow(isolatingSRGB)

% Get cone contrast image
for cc = 1:3
    temp = isolatingLMSImage(:,:,cc);
    meanLMSImage(cc) = mean(temp(:));
    temp = isolatingXYZImage(:,:,cc);
    meanXYZImage(cc) = mean(temp(:));
end
for ii = 1:imageN
    for jj = 1:imageN
        for cc = 1:3
            contrastLMSImage(ii,jj,cc) = (isolatingLMSImage(ii,jj,cc)-meanLMSImage(cc))/meanLMSImage(cc);
        end
    end
end
for cc = 1:3
    contrastLMS(cc,:) = (isolatingLMS(cc,:)-meanLMSImage(cc))/meanLMSImage(cc);
end

% Slice through LMS contrast image
% figure; hold on
% plot(1:imageN,gaborImage(centerN,:),'r','MarkerFaceColor','r+','MarkerSize',4);
% %plot(1:imageN,gaborImage(centerN+1,:),'g','MarkerFaceColor','g','MarkerSize',4);
% title('Gabor Image Slice');
% xlabel('x position (pixels)')
% ylabel('LMS Cone Contrast');
figure; hold on
plot(1:imageN,100*contrastLMSImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*contrastLMSImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*contrastLMSImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
title('Image Slice, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
% figure; hold on
% plot(1:imageN,isolatingImage(centerN,:,1),'ro','MarkerFaceColor','r','MarkerSize',4);
% plot(1:imageN,isolatingImage(centerN,:,2),'go','MarkerFaceColor','g','MarkerSize',4);
% plot(1:imageN,isolatingImage(centerN,:,3),'bo','MarkerFaceColor','b','MarkerSize',4);
% title('Image Slice, Rendered RGB');
% xlabel('x position (pixels)')
% ylabel('Rendered R,G,B');
% figure; hold on
% plot(1:imageN,isolatingXYZImage(centerN,:,1),'ro','MarkerFaceColor','r','MarkerSize',4);
% plot(1:imageN,isolatingXYZImage(centerN,:,2),'go','MarkerFaceColor','g','MarkerSize',4);
% plot(1:imageN,isolatingXYZImage(centerN,:,3),'bo','MarkerFaceColor','b','MarkerSize',4);
% title('Image Slice, Predicted XYZ');
% xlabel('x position (pixels)')
% ylabel('Rendered X, Y, Z');

%% Render colors
% theBgDeviceXYZ = T_xyz*theBgDeviceSpd;
% theIsolatingDeviceXYZUpper = T_xyz*theIsolatingDeviceSpdUpper;
% theIsolatingDeviceXYZLower = T_xyz*theIsolatingDeviceSpdLower;
% thBGDeviceSRGBPrimary = XYZToSRGBPrimary(theBgDeviceXYZ);
% theIsolatingDeviceSRGBPrimaryUpper = XYZToSRGBPrimary(theIsolatingDeviceXYZUpper);
% theIsolatingDeviceSRGBPrimaryLower  = XYZToSRGBPrimary(theIsolatingDeviceXYZLower );
% 
% scaleFactor = max([thBGDeviceSRGBPrimary ; theIsolatingDeviceSRGBPrimaryUpper; ; theIsolatingDeviceSRGBPrimaryLower]);
% theBGDeviceSRGB = SRGBGammaCorrect(thBGDeviceSRGBPrimary/scaleFactor,0);
% theIsolatingDeviceSRGBUpper = SRGBGammaCorrect(theIsolatingDeviceSRGBPrimaryUpper/scaleFactor,0);
% theIsolatingDeviceSRGBLower = SRGBGammaCorrect(theIsolatingDeviceSRGBPrimaryLower/scaleFactor,0);
% 
% theImage = zeros(imageN,imageN,3);
% for kk = 1:3
% 	theImage(:,:,kk) = theBGDeviceSRGB(kk);
% end
% for kk = 1:3
%     theImage((imageN-centerN)/2:imageN-(imageN-centerN)/2, ...
%              (imageN-centerN)/2:imageN-(imageN-centerN)/2, ...
%              kk) = theIsolatingDeviceSRGBUpper(kk);
% end
% figure; imshow(uint8(theImage));
% 
% theImage = zeros(imageN,imageN,3);
% for kk = 1:3
% 	theImage(:,:,kk) = theBGDeviceSRGB(kk);
% end
% for kk = 1:3
%     theImage((imageN-centerN)/2:imageN-(imageN-centerN)/2, ...
%              (imageN-centerN)/2:imageN-(imageN-centerN)/2, ...
%              kk) = theIsolatingDeviceSRGBLower(kk);
% end
% figure; imshow(uint8(theImage));

%% Some light level tests

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
theBGDeviceRawLum = T_xyz(2,:)*theBgDeviceSpd;
theBgDeviceSpdScaled = targetLum*theBgDeviceSpd/theBGDeviceRawLum;




