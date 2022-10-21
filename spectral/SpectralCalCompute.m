% SpectralCalCompute
%
% Explore spectral fits with swubprimaries, this
% version using the calibration structures.
%

% History:
%    04/22/2020  Started on it.

%% Clear
clear; close all;

%% Verbose?
%
% Set to true to get more output
VERBOSE = false;

% Set wavelength support.
%
% This needs to match what's in the calibration files, but
% we need it before we read those files.  A mismatch will
% throw an error below.
S = [380 2 201];

%% Set key stimulus parameters
%
% Condition Name.
conditionName = 'LminusMSmooth';
switch (conditionName)
    case 'LminusMSmooth'
        % Background xy.
        %
        % Specify the chromaticity, but we'll chose the luminance based
        % on the range available in the device.
        targetBgxy = [0.3127 0.3290]';

        % Target color direction and max contrasts.
        %
        % This is the basic desired modulation direction positive excursion. We go
        % equally in positive and negative directions.  Make this unit vector
        % length, as that is good convention for contrast.
        targetStimulusContrastDir = [1 -1 0]'; targetStimulusContrastDir = targetStimulusContrastDir/norm(targetStimulusContrastDir);

        % Specify desired primary properties.
        %
        % These are the target contrasts for the three primaries. We want these to
        % span a triangle around the line specified above. Here we define that
        % triangle by hand.  May need a little fussing for other directions, and
        % might be able to autocompute good choices.
        targetScreenPrimaryContrastDir(:,1) = [-1 1 0]'; targetScreenPrimaryContrastDir(:,1) = targetScreenPrimaryContrastDir(:,1)/norm(targetScreenPrimaryContrastDir(:,1));
        targetScreenPrimaryContrastDir(:,2) = [1 -1 0.5]'; targetScreenPrimaryContrastDir(:,2) = targetScreenPrimaryContrastDir(:,2)/norm(targetScreenPrimaryContrastDir(:,2));
        targetScreenPrimaryContrastDir(:,3) = [1 -1 -0.5]'; targetScreenPrimaryContrastDir(:,3) = targetScreenPrimaryContrastDir(:,3)/norm(targetScreenPrimaryContrastDir(:,3));

        % Set parameters for getting desired target primaries.
        targetScreenPrimaryContrasts = [0.05 0.05 0.05];
        targetPrimaryHeadroom = 1.05;
        primaryHeadroom = 0;
        targetLambda = 3;

        % We may not need the whole direction contrast excursion. Specify max
        % contrast we want relative to that direction vector.
        % The first number is
        % the amount we want to use, the second has a little headroom so we don't
        % run into numerical error at the edges. The second number is used when
        % defining the three primaries, the first when computing desired weights on
        % the primaries.
        spatialGaborTargetContrast = 0.04;
        plotAxisLimit = 100*spatialGaborTargetContrast;

        % Set up basis to try to keep spectra close to.
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
                B_naturalRaw = SplineSpd(S_cieday,B_cieday,S);
            case 'fourier'
                B_naturalRaw = MakeFourierBasis(S,nFourierBases);
            otherwise
                error('Unknown basis set specified');
        end
        B_natural{1} = B_naturalRaw;
        B_natural{2} = B_naturalRaw;
        B_natural{3} = B_naturalRaw;

    case 'ConeIsolating'
        % Background xy.
        %
        % Specify the chromaticity, but we'll chose the luminance based
        % on the range available in the device.
        targetBgxy = [0.3127 0.3290]';

        % Target color direction and max contrasts.
        %
        % This is the basic desired modulation direction positive excursion. We go
        % equally in positive and negative directions.  Make this unit vector
        % length, as that is good convention for contrast.
        targetStimulusContrastDir = [1 -1 0]'; targetStimulusContrastDir = targetStimulusContrastDir/norm(targetStimulusContrastDir);

        % Specify desired primary properties.
        %
        % These are the target contrasts for the three primaries. We want these to
        % span a triangle around the line specified above. Here we define that
        % triangle by hand.  May need a little fussing for other directions, and
        % might be able to autocompute good choices.
        targetScreenPrimaryContrastDir(:,1) = [1 0 0]'; targetScreenPrimaryContrastDir(:,1) = targetScreenPrimaryContrastDir(:,1)/norm(targetScreenPrimaryContrastDir(:,1));
        targetScreenPrimaryContrastDir(:,2) = [0 1 0]'; targetScreenPrimaryContrastDir(:,2) = targetScreenPrimaryContrastDir(:,2)/norm(targetScreenPrimaryContrastDir(:,2));
        targetScreenPrimaryContrastDir(:,3) = [0 0 1]'; targetScreenPrimaryContrastDir(:,3) = targetScreenPrimaryContrastDir(:,3)/norm(targetScreenPrimaryContrastDir(:,3));

        % Set parameters for getting desired target primaries.
        targetScreenPrimaryContrasts = [0.03 0.03 0.03];
        targetPrimaryHeadroom = 1;
        primaryHeadroom = 0.0;
        targetLambda = 3;

        % We may not need the whole direction contrast excursion. Specify max
        % contrast we want relative to that direction vector.
        % The first number is
        % the amount we want to use, the second has a little headroom so we don't
        % run into numerical error at the edges. The second number is used when
        % defining the three primaries, the first when computing desired weights on
        % the primaries.
        spatialGaborTargetContrast = 0.04;
        plotAxisLimit = 100*spatialGaborTargetContrast;

        % Set up basis to try to keep spectra close to.
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
                B_naturalRaw = SplineSpd(S_cieday,B_cieday,S);
            case 'fourier'
                B_naturalRaw = MakeFourierBasis(S,nFourierBases);
            otherwise
                error('Unknown basis set specified');
        end
        B_natural{1} = B_naturalRaw;
        B_natural{2} = B_naturalRaw;
        B_natural{3} = B_naturalRaw;
end

%% Define calibration filenames/params.
%
% This is a standard calibration file for the DLP screen,
% with the subprimaries set to something.  As we'll see below,
% we're going to rewrite those.nPrimaries
screenCalName = 'SACC';
screenNInputLevels = 256;

% These are the calibration files for each of the primaries, which
% then entails measuring the spectra of all the subprimaries for that
% primary.
channelCalNames = {'SACCPrimary1' 'SACCPrimary2' 'SACCPrimary3'};
channelNInputLevels = 253;

%% Load screen calibration and refit its gamma
screenCal = LoadCalFile(screenCalName);
screenCalObj = ObjectToHandleCalOrCalStruct(screenCal);
gammaMethod = 'identity';
screenCalObj.set('gamma.fitType',gammaMethod);
CalibrateFitGamma(screenCalObj, screenNInputLevels);

%% Load channel calibrations.
nScreenPrimaries = 3;
channelCals = cell(nScreenPrimaries ,1);
channelCalObjs = cell(nScreenPrimaries ,1);
for cc = 1:length(channelCalNames)
    channelCals{cc} = LoadCalFile(channelCalNames{cc});
    channelCalObjs{cc} = ObjectToHandleCalOrCalStruct(channelCals{cc});
    CalibrateFitGamma(channelCalObjs{cc}, channelNInputLevels);
end

%% Get out some data to work with.
%
% This is from the channel calibration file.
Scheck = channelCalObjs{1}.get('S');
if (any(S ~= Scheck))
    error('Mismatch between calibration file S and that specified at top');
end
wls = SToWls(S);
nChannels = channelCalObjs{1}.get('nDevices');

%% Cone fundamentals and XYZ CMFs.
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
psiParamsStruct.coneParams.fieldSizeDegrees = 2;
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
    gammaMeasurements = channelCals{1}.rawData.gammaCurveMeanMeasurements;
    [~,nMeas,~] = size(gammaMeasurements);
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

%% Image spatial parameters.
sineFreqCyclesPerImage = 6;
gaborSdImageFraction = 0.1;

% Image size in pixels
imageN = 512;

%% Get half on spectrum.
%
% This is useful for scaling things reasonably - we start with half of the
% available range of the primaries.
halfOnChannels = 0.5*ones(nChannels,1);
halfOnSpd = PrimaryToSpd(channelCalObjs{1},halfOnChannels);

%% Make sure gamma correction behaves well with unquantized conversion.
%
% This is just a check, not a computational step we need.
SetGammaMethod(channelCalObjs{1},0);
halfOnSettings = PrimaryToSettings(channelCalObjs{1},halfOnChannels);
halfOnPrimariesChk = SettingsToPrimary(channelCalObjs{1},halfOnSettings);
if (max(abs(halfOnChannels-halfOnPrimariesChk)) > 1e-8)
    error('Gamma self-inversion not sufficiently precise');
end

%% Use quantized conversion from here on.
%
% Comment in the line that refits the gamma to see
% effects of extreme quantization one what follows.
%
% CalibrateFitGamma(channelCalObjs{1},10);
channelGammaMethod = 2;
SetGammaMethod(channelCalObjs{1},channelGammaMethod);
SetGammaMethod(channelCalObjs{2},channelGammaMethod);
SetGammaMethod(channelCalObjs{3},channelGammaMethod);

%% Use extant machinery to get primaries from spectrum.
%
% This isn't used in our calculations.  Any difference in the
% two lines here reflects a bug in the SpdToPrimary/PrimaryToSpd pair.
if (VERBOSE)
    halfOnPrimariesChk = SpdToPrimary(channelCalObjs{1},halfOnSpd);
    halfOnSpdChk = PrimaryToSpd(channelCalObjs{1},halfOnPrimariesChk);
    figure; hold on;
    plot(wls,halfOnSpd,'r','LineWidth',3);
    plot(wls,halfOnSpdChk,'k','LineWidth',1);

    % Show effect of quantization.
    %
    % It's very small at the nominal 253 levels of the subprimaries, but will
    % increase if you refit the gamma functios to a small number of levels.
    halfOnPrimariesChk = SpdToPrimary(channelCalObjs{1},halfOnSpd);
    halfOnSettingsChk = PrimaryToSettings(channelCalObjs{1},halfOnPrimariesChk);
    halfOnPrimariesChk1 = SettingsToPrimary(channelCalObjs{1},halfOnSettingsChk);
    halfOnSpdChk1 = PrimaryToSpd(channelCalObjs{1},halfOnPrimariesChk1);
    plot(wls,halfOnSpdChk1,'g','LineWidth',1);
end

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
targetBgXYZ = xyYToXYZ([targetBgxy ; 1]);

% Adjust these to keep background in gamut
primaryBackgroundScaleFactor = 0.5;
screenBackgroundScaleFactor = 0.5;

% Make a loop for getting background for all primaries.
% Passing true for key 'Scale' causes these to be scaled reasonably
% relative to gamut, which is why we can set the target luminance
% arbitrarily to 1 just above. The scale factor determines where in the
% approximate channel gamut we aim the background at.
for pp = 1:nScreenPrimaries
    [channelBackgroundPrimaries(:,pp),channelBackgroundSpd(:,pp),channelBackgroundXYZ(:,pp)] = FindBgChannelPrimaries(targetBgXYZ,T_xyz,channelCalObjs{pp}, ...
        B_natural{pp},projectIndices,primaryHeadRoom,targetLambda,'scaleFactor',0.6,'Scale',true,'Verbose',true);
end
if (any(channelBackgroundPrimaries < 0) | any(channelBackgroundPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end
fprintf('Background primary min: %0.2f, max: %0.2f, mean: %0.2f\n', ...
    min(channelBackgroundPrimaries(:)),max(channelBackgroundPrimaries(:)),mean(channelBackgroundPrimaries(:)));

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
        = FindChannelPrimaries(targetScreenPrimaryContrastDir(:,pp), ...
        targetPrimaryHeadroom,targetScreenPrimaryContrasts(pp),channelBackgroundPrimaries(:,pp), ...
        T_cones,channelCalObjs{pp},B_natural{pp},projectIndices,primaryHeadroom,targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
    
    % We can wonder about how close to gamut our primaries are.  Compute
    % that here.
    primaryGamutScaleFactor(pp) = MaximizeGamutContrast(screenPrimaryModulationPrimaries(:,pp),channelBackgroundPrimaries(:,pp));
    fprintf('\tPrimary %d, gamut scale factor is %0.3f\n',pp,primaryGamutScaleFactor(pp));
    
    % Find the channel settings that correspond to the desired screen
    % primaries.
    screenPrimarySettings(:,pp) = PrimaryToSettings(channelCalObjs{pp},screenPrimaryPrimaries(:,pp));
end

%% How close are spectra to subspace defined by basis?
isolatingNaturalApproxSpd1 = B_natural{1}*(B_natural{1}(projectIndices,:)\screenPrimarySpd(projectIndices,1));
isolatingNaturalApproxSpd2 = B_natural{2}*(B_natural{2}(projectIndices,:)\screenPrimarySpd(projectIndices,2));
isolatingNaturalApproxSpd3 = B_natural{3}*(B_natural{3}(projectIndices,:)\screenPrimarySpd(projectIndices,3));

% Plot of the screen primary spectra.
subplot(2,2,1); hold on
plot(wls,screenPrimarySpd(:,1),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd1,'r:','LineWidth',1);
plot(wls(projectIndices),screenPrimarySpd(projectIndices,1),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd1(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 1');

subplot(2,2,2); hold on
plot(wls,screenPrimarySpd(:,2),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd2,'r:','LineWidth',1);
plot(wls(projectIndices),screenPrimarySpd(projectIndices,2),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd2(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 2');

subplot(2,2,3); hold on
plot(wls,screenPrimarySpd(:,3),'b','LineWidth',2);
plot(wls,isolatingNaturalApproxSpd3,'r:','LineWidth',1);
plot(wls(projectIndices),screenPrimarySpd(projectIndices,3),'b','LineWidth',4);
plot(wls(projectIndices),isolatingNaturalApproxSpd3(projectIndices),'r:','LineWidth',3);
xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
title('Primary 3');

%% Set the screen primaries.
%
% We want these to match those we set up with the
% channel calculations above.  Need to reset
% sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device',screenPrimarySpd);
SetSensorColorSpace(screenCalObj,T_cones,S);

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

%% Set up desired background.
%
% We aim for the background that we said we wanted when we built the screen primaries.
desiredBgExcitations = screenBackgroundScaleFactor*T_cones*sum(channelBackgroundSpd,2);
screenBgSettings = SensorToSettings(screenCalObj,desiredBgExcitations);
screenBgExcitations = SettingsToSensor(screenCalObj,screenBgSettings);
figure; clf; hold on;
plot(desiredBgExcitations,screenBgExcitations,'ro','MarkerFaceColor','r','MarkerSize',12);
axis('square');
xlim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
ylim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
xlabel('Desired bg excitations'); ylabel('Obtained bg excitations');
title('Check that we obtrain desired background excitations');
fprintf('Screen settings to obtain background: %0.2f, %0.2f, %0.2f\n', ...
    screenBgSettings(1),screenBgSettings(2),screenBgSettings(3));

%% Make monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image.  Multiply
% by the max contrast vector to get the LMS contrast image.
fprintf('Making Gabor contrast image\n');
centerN = imageN/2;
gaborSdPixels = gaborSdImageFraction*imageN;
rawMonochromeSineImage = MakeSineImage(0,sineFreqCyclesPerImage,imageN);
gaussianWindow = normpdf(MakeRadiusMat(imageN,imageN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));
rawMonochromeUnquantizedContrastGaborImage = rawMonochromeSineImage.*gaussianWindow;

% Put it into cal format.  Each pixel in cal format is one column.  Here
% there is just one row since it is a monochrome image at this point.
rawMonochromeUnquantizedContrastGaborCal = ImageToCalFormat(rawMonochromeUnquantizedContrastGaborImage);

%% Quantize the contrast image to a (large) fixed number of levels.
%
% This allows us to speed up the image conversion without any meaningful
% loss of precision. If you don't like it, increase number of quantization
% bits until you are happy again.
nQuantizeBits = 14;
nQuantizeLevels = 2^nQuantizeBits;
rawMonochromeContrastGaborCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastGaborCal+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;

% Plot of how well point cloud method does in obtaining desired contrats.
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
desiredContrastGaborCal = spatialGaborTargetContrast*targetStimulusContrastDir*rawMonochromeContrastGaborCal;

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

%% Find settings by exhaustive search of point cloud for each pixel
%
% Go through the gabor image, and for each pixel find the settings that
% come as close as possible to producing the desired excitations.
% Conceptually straightforward, but a bit slow.
SLOWMETHODCHECK = false;
if (SLOWMETHODCHECK)
    tic;
    fprintf('Point cloud exhaustive method, finding image settings\n')
    printIter = 10000;
    ptCldSettingsGaborCal = zeros(3,size(desiredContrastGaborCal,2));
    minIndex = zeros(1,size(desiredContrastGaborCal,2));
    for ll = 1:size(desiredContrastGaborCal,2)
        if (rem(ll,printIter) == 0)
            fprintf('Finding settings for iteration %d of %d\n',ll,size(desiredContrastGaborCal,2));
        end
        minIndex = findNearestNeighbors(contrastPtCld,desiredContrastGaborCal(:,ll)',1);
        ptCldSettingsGaborCal(:,ll) = ptCldSettingsCal(:,minIndex);
    end
    toc

    % Get contrasts we think we have obtained.
    ptCldExcitationsGaborCal = SettingsToSensor(screenCalObj,ptCldSettingsGaborCal);
    ptCldContrastGaborCal = ExcitationsToContrast(ptCldExcitationsGaborCal,screenBgExcitations);

    % Plot of how well pixelwise point cloud method does in obtaining desired contrats
    figure; clf;
    plot(desiredContrastGaborCal(:),ptCldContrastGaborCal(:),'r+');
    axis('square');
    xlabel('Desired L, M or S contrast');
    ylabel('Predicted L, M, or S contrast');
    title('Pixelwise point cloud image method');
end

%% Get image from point cloud, in cal format
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

% Check that we get the same answer
if (SLOWMETHODCHECK)
    if (max(abs(uniqueQuantizedContrastGaborCal(:)-uniqueQuantizedContrastGaborCal(:))) > 0)
        fprintf('Point cloud and unique method methods do not agree\n');
    end
end

%% Convert representations we want to take forward to image format
desiredContrastGaborImage = CalFormatToImage(desiredContrastGaborCal,imageN,imageN);
standardPredictedContrastImage = CalFormatToImage(standardPredictedContrastGaborCal,imageN,imageN);
standardSettingsGaborImage = CalFormatToImage(standardSettingsGaborCal,imageN,imageN);
uniqueQuantizedContrastGaborImage = CalFormatToImage(uniqueQuantizedContrastGaborCal,imageN,imageN);

%% SRGB image via XYZ, scaled to display
predictedXYZCal = T_xyz*desiredSpdGaborCal;
SRGBPrimaryCal = XYZToSRGBPrimary(predictedXYZCal);
scaleFactor = max(SRGBPrimaryCal(:));
SRGBCal = SRGBGammaCorrect(SRGBPrimaryCal/(2*scaleFactor),0);
SRGBImage = uint8(CalFormatToImage(SRGBCal,imageN,imageN));

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
figure; hold on
plot(1:imageN,100*standardPredictedContrastImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*standardPredictedContrastImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*standardPredictedContrastImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
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
plot(1:imageN,100*uniqueQuantizedContrastGaborImage(centerN,:,1),'r+','MarkerFaceColor','r','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:imageN,100*uniqueQuantizedContrastGaborImage(centerN,:,2),'g+','MarkerFaceColor','g','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:imageN,100*uniqueQuantizedContrastGaborImage(centerN,:,3),'b+','MarkerFaceColor','b','MarkerSize',4);
plot(1:imageN,100*desiredContrastGaborImage(centerN,:,3),'b','LineWidth',0.5);
title('Image Slice, Point Cloud Method, LMS Cone Contrast');
xlabel('x position (pixels)')
ylabel('LMS Cone Contrast (%)');
ylim([-plotAxisLimit plotAxisLimit]);

%% Generate some settings values corresponding to known contrasts
%
% The reason for this is to measure and check these.  This logic follows
% how we handled an actual gabor image above. We don't actually need to
% quantize to 14 bits here on the contrast, but nor does it hurt.
rawMonochromeUnquantizedContrastCheckCal = [0 0.05 -0.05 0.10 -0.10 0.15 -0.15 0.20 -0.20 0.25 -0.25 0.5 -0.5 1 -1];
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
desiredContrastCheckCal = spatialGaborTargetContrast*targetStimulusContrastDir*rawMonochromeContrastCheckCal;
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
plot(SToWls(S),screenPrimarySpdCheck,'k','LineWidth',4);
plot(SToWls(S),screenPrimarySpd,'r','LineWidth',2);
xlabel('Wavelength'); ylabel('Radiance');
title('Check of consistency between screen primaries and screen primary spds');

%% Save out what we need to check things on the DLP
screenSettingsImage = standardSettingsGaborImage;
if (ispref('SpatioSpectralStimulator','SACCData'))
    dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    testFiledir = getpref('SpatioSpectralStimulator','SACCData');
    testFilename = fullfile(testFiledir,'CheckCalibration',sprintf('testImageData_%s',dayTimestr));
    save(testFilename,'S','T_cones','screenCalObj','channelCalObjs','screenSettingsImage', ...
        'screenPrimaryPrimaries','screenPrimarySettings','screenPrimarySpd',...
        'desiredContrastCheckCal', ...
        'ptCldScreenSettingsCheckCal','ptCldScreenContrastCheckCal','ptCldScreenSpdCheckCal', ...
        'nQuantizeLevels','screenNInputLevels','targetStimulusContrastDir','spatialGaborTargetContrast');
end
