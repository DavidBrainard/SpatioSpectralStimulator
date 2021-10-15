% SpectralTestCal
%
% Start exploring spectral fits with swubprimarys, this
% version using the calibration structures.
%
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Verbose
%
% Set to true to get more output
VERBOSE = false;

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
% equally in positive and negative directions.  Make this unit vector
% length, as that is good convention for contrast.
targetContrastDir = [1 -1 0]';
targetContrastDir = targetContrastDir/norm(targetContrastDir);

% We may not need the whole direction contrast excursion. Specify max
% contrast we want relative to that direction vector.
% The first number is
% the amount we want to use, the second has a little headroom so we don't
% run into numerical error at the edges. The second number is used when
% defining the three primaries, the first when computing desired weights on
% the primaries.
targetContrast = 0.05;
plotAxisLimit = 100*targetContrast;

%% Specify desired primary properties.
%
% These are the target contrasts for the three primaries. We want these to
% span a triangle around the line specified above. Here we define that
% triangle by hand.  May need a little fussing for other directions, and
% might be able to autocompute good choices.
targetPrimary1MaxContrast = [-1 1 0]'; targetPrimary1MaxContrast = targetPrimary1MaxContrast/norm(targetPrimary1MaxContrast);
targetPrimary2MaxContrast = [1 -1 0.5]'; targetPrimary2MaxContrast = targetPrimary2MaxContrast/norm(targetPrimary2MaxContrast);
targetPrimary3MaxContrast = [1 -1 -0.5]'; targetPrimary3MaxContrast = targetPrimary3MaxContrast/norm(targetPrimary3MaxContrast);

%% Image spatial parameters.
sineFreqCyclesPerImage = 6;
gaborSdImageFraction = 0.1;

% Image size in pixels
imageN = 512;

%% Computational bit depth.
%
% This is a computational bit depth that we use to define the lookup table
% between contrast and primary values.
% fineBits = 14;
% nFineLevels = 2^fineBits;

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
subprimaryGammaMethod = 2;
SetGammaMethod(subprimaryCalObjs{1},subprimaryGammaMethod);
SetGammaMethod(subprimaryCalObjs{2},subprimaryGammaMethod);
SetGammaMethod(subprimaryCalObjs{3},subprimaryGammaMethod);

%% Use extant machinery to get primaries from spectrum.
%
% This isn't used in our calculations.  Any difference in the
% two lines here reflects a bug in the SpdToPrimary/PrimaryToSpd pair.
if (VERBOSE)
    halfOnPrimariesChk = SpdToPrimary(subprimaryCalObjs{1},halfOnSpd);
    halfOnSpdChk = PrimaryToSpd(subprimaryCalObjs{1},halfOnPrimariesChk);
    figure; hold on;
    plot(wls,halfOnSpd,'r','LineWidth',3);
    plot(wls,halfOnSpdChk,'k','LineWidth',1);

    % Show effect of quantization.
    %
    % It's very small at the nominal 252 levels of the subprimaries, but will
    % increase if you refit the gamma functios to a small number of levels.
    halfOnPrimariesChk = SpdToPrimary(subprimaryCalObjs{1},halfOnSpd);
    halfOnSettingsChk = PrimaryToSettings(subprimaryCalObjs{1},halfOnPrimariesChk);
    halfOnPrimariesChk1 = SettingsToPrimary(subprimaryCalObjs{1},halfOnSettingsChk);
    halfOnSpdChk1 = PrimaryToSpd(subprimaryCalObjs{1},halfOnPrimariesChk1);
    plot(wls,halfOnSpdChk1,'g','LineWidth',1);
end

%% Set up basis to try to keep spectra close to.
%
% This is how we enforce a smoothness or other constraint
% on the spectra.  What happens in the routine that finds
% primaries is that there is a weighted error term that tries to
% maximize the projection onto a passed basis set.
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
%
% Set parameters for getting desired background primaries.
primaryHeadRoom = 0;
targetLambda = 3;
targetBgXYZ = xyYToXYZ([targetBgxy ; 1]);

% Make a loop for getting background for all primaries.
% Passing true for key 'Scale' causes these to be scaled reasonably
% relative to gamut, which is why we can set the target luminance
% arbitrarily to 1 just above.
for pp = 1:nPrimaries
    [bgPrimaries(:,pp),obtainedBgSpd(:,pp),obtainedBgXYZ(:,pp)] = FindDesiredBackgroundPrimaries(targetBgXYZ,T_xyz,subprimaryCalObjs{pp}, ...
        B_natural,projectIndices,primaryHeadRoom,targetLambda,'Scale',true,'Verbose',true);
end
if (any(bgPrimaries < 0) | any(bgPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end

%% Find primaries with desired LMS contrast.
%
% Set parameters for getting desired target primaries.
targetMaxLMSContrast = [targetPrimary1MaxContrast targetPrimary2MaxContrast targetPrimary3MaxContrast];
targetPrimaryContrast = 0.07;
targetPrimaryHeadroom = 1.05;
primaryHeadroom = 0;
targetLambda = 3;

% Make a loop for getting isolating primaries for all primaries.
for pp = 1:nPrimaries

    % The ambient with respect to which we compute contrast is from all
    % three primaries, which we handle via the extraAmbientSpd key-value
    % pair in the call.  The extra is for the primaries not being found in
    % the current call - the contribution from the current primary is known
    % because we pass the primaries for the background.
    otherPrimaries = setdiff(1:nPrimaries,pp);
    extraAmbientSpd = 0;
    for oo = 1:length(otherPrimaries)
        extraAmbientSpd = extraAmbientSpd + obtainedBgSpd(:,otherPrimaries(oo));
    end

    % Get isolating primaries.
    [isolatingPrimaries(:,pp),isolatingPrimariesQuantized(:,pp),isolatingSpd(pp,:),isolatingContrast(pp,:)] = FindDesiredContrastTargetPrimaries(targetMaxLMSContrast(:,pp), ...
        targetPrimaryHeadroom,targetPrimaryContrast,bgPrimaries(:,pp), ...
        T_cones,subprimaryCalObjs{pp},B_natural,projectIndices,primaryHeadroom,targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
end

%% Measure the desired target primaries
%
% This result will be used to update the calibration file
% so that it matches reality rather than aspiration.
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

%% Plot of the background and primary spectra
figure; clf;
subplot(2,2,1); hold on
plot(wls,obtainedBgSpd,'b','LineWidth',2);
plot(wls,theBgNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),obtainedBgSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theBgNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Background');

subplot(2,2,2); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(1,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd1,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(1,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd1(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 1');

subplot(2,2,3); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(2,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd2,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(2,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd2(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 2');

subplot(2,2,4); hold on
plot(wls,obtainedBgSpd,'b:','LineWidth',1);
plot(wls,isolatingSpd(3,:),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd3,'r:','LineWidth',1);
plot(wls(projectIndices),isolatingSpd(3,projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd3(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 3');

%% Set the projector primaries
%
% We want these to match those we set up with the
% subprimary calculations above.  Need to reset
% sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
MEASURED = false;
if (~MEASURED)
    projectorCalObj.set('P_device',isolatingSpd');
else
    projectorCalObj.set('P_device',isolatingSpdMeasured);
end
SetSensorColorSpace(projectorCalObj,T_cones,S);

%% Set projector gamma method
%
% If we set to 0, there is no quantization and the result is excellent.
% If we set to 2, this is quantized at 256 levels and the result is more
% of a mess.  The choice of 2 represents what we think will actually happen
% since the real device is quantized.
%
% The point cloud method below reduces this problem.
projectorGammaMethod = 2;
SetGammaMethod(projectorCalObj,projectorGammaMethod);

%% Make Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image.  Multiply
% by the max contrast vector to get the LMS contrast image.
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawMonochromeSineImage = MakeSineImage(0,sineFreqCyclesPerImage,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawMonochromeContrastGaborImage = rawMonochromeSineImage.*gaussianWindow;

% Put it into cal format.  Each pixel in cal format is one column.  Here
% there is just one row since it is a monochrome image at this point.
rawMonochromeContrastGaborCal = ImageToCalFormat(rawMonochromeContrastGaborImage);

% Convert background to LMS
%
% Background is sum of backgrounds for each of the three primaries
projectorBgSettings = [0.5 0.5 0.5]';
projectorBgExcitations = SettingsToSensor(projectorCalObj,[projectorBgSettings]);

% Scale target LMS contrast vector at max excursion by contrast modulation
% at each pixel.  This is done by a single matrix multiply plus a lead
% factor.  We work cal format here as that makes color transforms
% efficient.
theDesiredContrastGaborCal = targetContrast*targetContrastDir*rawMonochromeContrastGaborCal;

% Convert cone contrast to excitations
theDesiredExcitationsGaborCal = ContrastToExcitation(theDesiredContrastGaborCal,projectorBgExcitations);

% Get primaries using standard calibration code
thePrimariesGaborCal = SensorToPrimary(projectorCalObj,theDesiredExcitationsGaborCal);

% Gamma correct and quantize (if gamma method set to 2 above; with gamma
% method set to zero there is no quantization).  Then convert back from
% the gamma corrected settings.
theSettingsGaborCal = PrimaryToSettings(projectorCalObj,thePrimariesGaborCal);
thePredictedPrimariesGaborCal = SettingsToPrimary(projectorCalObj,theSettingsGaborCal);
thePredictedExcitationsGaborCal = PrimaryToSensor(projectorCalObj,thePredictedPrimariesGaborCal);
thePredictedContrastGaborCal = ExcitationsToContrast(thePredictedExcitationsGaborCal,projectorBgExcitations);
thePredictedSpdGaborCal = PrimaryToSpd(projectorCalObj,thePredictedPrimariesGaborCal);

%% Do SensorToSettings by finding nearest points in point cloud
%
% The method above is subject to imperfect quantization because each primary is
% quantized individually. Here we'll quantize jointly across the three
% primaries, using an exhaustive search process.  Amazingly, it is feasible
% to search all possible quantized settings for each image pixel, and choose
% the settings that best approximate the desired LMS excitations at that pixel.
%
% Compute an array with all possible triplets of projector settings,
% quantized on the interval [0,1].
%
% This method takes all possible projector settings and creates a
% point cloud of the corresponding cone contrasts. It then finds
% the settings that come as close as possible to producing the
% desired cone contrast at each point in the image. It's a little
% slow but conceptually simple and fast enough to be feasible.
tic;
fprintf('Point cloud exhaustive method, setting up cone contrast cloud, this takes a while\n')
allProjectorSettingsCal = zeros(3,projectorNInputLevels^3);
idx = 1;
for ii = 0:(projectorNInputLevels-1)
    for jj = 0:(projectorNInputLevels-1)
        for kk = 0:(projectorNInputLevels-1)
            allProjectorSettingsCal(:,idx) = [ii jj kk]'/projectorNInputLevels;
            idx = idx+1;
        end
    end
end

% Get LMS excitations for each triplet of projector settings, and build a
% point cloud object from these.
allProjectorExcitationsGaborCal = SettingsToSensor(projectorCalObj,allProjectorSettingsCal);
allProjectorContrastGaborCal = ExcitationsToContrast(allProjectorExcitationsGaborCal,projectorBgExcitations);
allSensorPtCloud = pointCloud(allProjectorContrastGaborCal');

% Force point cloud setup by finding one nearest neighbor. This is slow,
% but once it is done subsequent calls are considerably faster.
findNearestNeighbors(allSensorPtCloud,[0 0 0],1);
toc

% Go through the gabor image, and for each pixel find the settings that
% come as close as possible to producing the desired excitations.
tic;
fprintf('Point cloud exhaustive method, finding image settings\n')
printIter = 10000;
thePointCloudSettingsGaborCal = zeros(3,size(theDesiredContrastGaborCal,2));
minIndex = zeros(1,size(theDesiredContrastGaborCal,2));
for ll = 1:size(theDesiredContrastGaborCal,2)
    if (rem(ll,printIter) == 0)
        fprintf('Finding settings for iteration %d of %d\n',ll,size(theDesiredContrastGaborCal,2));
    end
    minIndex(ll) = findNearestNeighbors(allSensorPtCloud,theDesiredContrastGaborCal(:,ll)',1);
    thePointCloudSettingsGaborCal(:,ll) = allProjectorSettingsCal(:,minIndex(ll));
end

% Get contrasts we think we have obtianed.
thePointCloudExcitationsGaborCal = SettingsToSensor(projectorCalObj,thePointCloudSettingsGaborCal);
thePointCloudContrastGaborCal = ExcitationsToContrast(thePointCloudExcitationsGaborCal,projectorBgExcitations);
toc

% Plot of how well point cloud method does in obtaining desired contrats
figure; clf;
plot(theDesiredContrastGaborCal(:),thePointCloudContrastGaborCal(:),'r+');
axis('square');
xlabel('Desired L, M or S contrast');
ylabel('Predicted L, M, or S contrast');

%% Convert representations we want to take forward to image format
theDesiredLMSContrastGaborImage = CalFormatToImage(theDesiredContrastGaborCal,imageN,imageN);
thePredictedLMSContrastImage = CalFormatToImage(thePredictedContrastGaborCal,imageN,imageN);
theSettingsGaborImage = CalFormatToImage(theSettingsGaborCal,imageN,imageN);
thePointCloudLMSContrastGaborImage = CalFormatToImage(thePointCloudContrastGaborCal,imageN,imageN);

%% SRGB image via XYZ, scaled to display
thePredictedXYZCal = T_xyz*thePredictedSpdGaborCal;
theSRGBPrimaryCal = XYZToSRGBPrimary(thePredictedXYZCal);
scaleFactor = max(theSRGBPrimaryCal(:));
theSRGBCal = SRGBGammaCorrect(theSRGBPrimaryCal/(2*scaleFactor),0);
theSRGBImage = uint8(CalFormatToImage(theSRGBCal,imageN,imageN));

% Show the SRGB image
figure; imshow(theSRGBImage);

%% Show the settings image
figure; clf;
imshow(theSettingsGaborImage);

%% Plot slice through predicted LMS contrast image.
%
% Note that the y-axis in this plot is individual cone contrast, which is
% not the same as the vector length contrast of the modulation.
figure; hold on
plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*thePredictedLMSContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
if (projectorGammaMethod == 2)
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
plot(1:imageN,100*thePointCloudLMSContrastGaborImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*thePointCloudLMSContrastGaborImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*thePointCloudLMSContrastGaborImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*theDesiredLMSContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, Point Cloud Method, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Generate some settings values corresponding to known contrasts
%
% The reason for this is to measure and check these.  This logic follows
% how we handled an actual gabor image above.
rawMonochromeContrastCheckCal = [0 0.25 -0.25 0.5 -0.5 1 -1];
theDesiredContrastCheckCal = targetContrast*targetContrastDir*rawMonochromeContrastCheckCal;
theDesiredExcitationsCheckCal = ContrastToExcitation(theDesiredContrastCheckCal,projectorBgExcitations);

% For each pixel find the settings that
% come as close as possible to producing the desired excitations.
%
% If we measure for a uniform field the spectra corresopnding to each of
% the settings in the columns of thePointCloudSettingsCheckCal, then
% compute the cone contrasts with respect to the backgound (0 contrast
% measurement, first settings), we should approximate the cone contrasts in
% theDesiredContrastCheckCal.
fprintf('Point cloud exhaustive method, finding check settings\n')
thePointCloudSettingsCheckCal = zeros(3,size(theDesiredContrastCheckCal,2));
for ll = 1:size(theDesiredContrastCheckCal,2)
    minIndex = findNearestNeighbors(allSensorPtCloud,theDesiredContrastCheckCal(:,ll)',1);
    thePointCloudSettingsCheckCal(:,ll) = allProjectorSettingsCal(:,minIndex);
end
thePointCloudExcitationsCheckCal = SettingsToSensor(projectorCalObj,thePointCloudSettingsCheckCal);
thePointCloudContrastCheckCal = ExcitationsToContrast(thePointCloudExcitationsCheckCal,projectorBgExcitations);

%% Save out what we need to check things on the DLP
projectorSettingsImage = theSettingsGaborImage;
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,'testImageData1');
    save(testFilename,'projectorSettingsImage','isolatingPrimaries','theDesiredContrastCheckCal','thePointCloudSettingsCheckCal','thePointCloudContrastCheckCal');
end



