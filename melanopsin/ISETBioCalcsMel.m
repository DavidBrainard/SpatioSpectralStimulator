%{

% Asano et al. give the following population SD's for the individual
% difference parameters (their Table 5, Step 2 numbers:
%   Lens    - 18.7%
%   Macular - 36.5%
%   L Density - 9%
%   M Density - 9%
%   S Density - 7.4%
%   L Shift   - 2 nm
%   M Shift   - 1.5 nm
%   S Shift   - 1.3 nm

%% Run out a bunch of examples
clear;
conditionNameList = {'MelDirected1' 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2 1 2 5 10];
gaborSdDeg = 100;
stimulusSizeDeg = 4;

fieldSizeDeg = 5;
eccXDegList = [0 -5 -10 -15 -20];
aoRender = false;
noLCA = false;
ageISETBioList = [20 32 60];
dLensISETBioList = [-18.7 0 18.7];
dMacISETBioList = [-36.5 0 36.5];

for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        for ee = 1:length(eccXDegList)
            for aa = 1:length(ageISETBioList)
                for ll = 1:length(dLensISETBioList)
                    for mm = 1:length(dMacISETBioList)
                        lmsConeContrast(:,cc,ss,ee,aa,ll,mm) = ISETBioCalcsMel(conditionNameList{cc},sineFreqCyclesPerDegList(ss), ...
                                gaborSdDeg,stimulusSizeDeg, ...
                                fieldSizeDeg, ...
                                eccXDegList(ee), ...
                                aoRender, noLCA, ...
                                ageISETBioList(aa), dLensISETBioList(ll), dMacISETBioList(mm) ...
                            );
                    end
                end
            end
        end
    end
end

%}

%{

%% Baseline condition test
clear;
conditionName = 'MelDirected1';
sineFreqCyclesPerDeg = 0.2;
gaborSdDeg = 100;
stimulusSizeDeg = 4;

fieldSizeDeg = 5;
eccXDeg = -5;
aoRender = false;
noLCA = false;
ageISETBio = 60;
dLensISETBio = -5;
dMacISETBio = 30;

[lmsConeContrast] = ISETBioCalcsMel(conditionName,sineFreqCyclesPerDeg, ...
        gaborSdDeg,stimulusSizeDeg, ...
        fieldSizeDeg, ...
        eccXDeg, ...
        aoRender, noLCA, ...
        ageISETBio, dLensISETBio, dMacISETBio ...
    );

%}

% ISETBioCalcsMel
function [lmsConeContrast] = ISETBioCalcsMel(conditionName,sineFreqCyclesPerDeg, ...
            gaborSdDeg,stimulusSizeDeg, ...
            fieldSizeDeg, ...
            eccXDeg, ...
            aoRender,noLCA, ...
            ageISETBio, dLensISETBio, dMacISETBio ...
        )

%% Initialize
close all;

%% Parameter specifications for scene that are not passed.
screenGammaMethod = 2;
sceneInputStr = sprintf('%s_Size_%0.1f_Sf_%0.1f_Sd_%0.1f_GammaMethod_%d', ...
    conditionName,stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenGammaMethod);

% Load the scene data according to parameters above
projectFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
sceneInputSubdir = fullfile(projectFiledir,sceneInputStr);
sceneInputFilename = fullfile(sceneInputSubdir,'sceneOutput.mat');
theData = load(sceneInputFilename);
disp('Data loaded');

%% Mosaic settings not so likely to vary
eccYDeg = 0;
subjectID = 6;
mosaicOutputStr = sprintf('%s_EccX_%0.1f_EccY_%0.1f_FieldSize_%0.1f_AO_%d_NoLCA_%d_Age_%d_macAdj_%0.1f_lensAdj_%0.1f', ...
    conditionName,eccXDeg,eccYDeg,fieldSizeDeg,aoRender,noLCA,ageISETBio,dMacISETBio,dLensISETBio);

%% Use ISETPipelineToolbox wrapper as a way to the isetbio computations
%
% These parameters define key addtional properties we are likely to vary
% Not likely to vary these in the short run, so not coded in output name for now
randSeed = false;                              % False means, no randomness
eccVars = false;
pupilDiamMM = 3;
defocusDiopters = 0;
zernikeDataBase = 'Polans2015';

%% Setup output names and places
outputDir = fullfile(sceneInputSubdir,mosaicOutputStr);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
mosaicOutputFilename = fullfile(outputDir,[sprintf('mosaicOutput') '.mat']);

%% Compute radial eccentricity
eccDeg = sqrt(eccXDeg^2 + eccYDeg^2);

%% Generate the photopigment, macular pigment density, and lens density to use
%
% What we're going to do is override the oi and cMosaic's view of what
% these are, so we can control them completely.  We'll calculate these
% first so that we can insert the cMosaic parameters on creation, which is
% easier than doing it ex post.  See tutorial t_conesPhotopigment.
%
% Lens to use
lensTransmittance = LensTransmittance(theData.colorDirectionParams.wls,'Human','CIE', ageISETBio, pupilDiamMM);
unadjustedLensDensity = -log10(lensTransmittance);
adjustedLensDensity = unadjustedLensDensity * (1 + dLensISETBio/100);
lensTransmittance = 10.^-adjustedLensDensity;
lensObject = Lens('wave',theData.colorDirectionParams.wls,'unitDensity',-log10(lensTransmittance),'density',1);

% Plot of lens density versus age
plotAges = linspace(20,80,100);
for aa = 1:length(plotAges)
    lensTransmittanceAA = LensTransmittance(theData.colorDirectionParams.wls,'Human','CIE', plotAges(aa), pupilDiamMM);
    unadjustedLensDensityAA = -log10(lensTransmittanceAA);
    plotLensDensity(aa) = max(unadjustedLensDensityAA(:));
end
figure; clf; hold on
plot(plotAges,plotLensDensity,'b','LineWidth',6);
xlabel('Age (years)','FontName','Helvetica','FontSize',18);
ylabel('Lens Peak Density','FontName','Helvetica','FontSize',18);
ylim([0 6]);
saveas(gcf,fullfile(outputDir,'lensDensityWithAge.tiff'),'tiff');

%% Macular transmittance for the specified field size
macTransmittance = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',2*eccDeg);
unadjustedMacDensity = -log10(macTransmittance);
adjustedMacDensity = unadjustedMacDensity * (1 + dMacISETBio/100);
macTransmittance = 10.^-adjustedMacDensity;
macObject = Macular('wave',theData.colorDirectionParams.wls,'unitDensity',-log10(macTransmittance),'density',1);

% Plot of macular pigment density versus eccentricity
plotEccs = linspace(0,25,100);
for ee = 1:length(plotEccs)
    macTransmittanceEE = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',2*plotEccs(ee));
    unadjustedMacDensityEE = -log10(macTransmittanceEE);
    plotMacDensity(ee) = max(unadjustedMacDensityEE(:));
end
figure; clf; hold on
plot(plotEccs,plotMacDensity,'b','LineWidth',6);
xlabel('Eccentricity (degs)','FontName','Helvetica','FontSize',18);
ylabel('Macular Pigment Peak Density','FontName','Helvetica','FontSize',18);
ylim([0 max(plotMacDensity(:))]);
saveas(gcf,fullfile(outputDir,'macDensityWithEccentricty.tiff'),'tiff');

%% Photopigment absorptance
%
% The individual differences model allows for shifting the absorbance along
% the wavelength axis. Do that here.
asanoConeParams = DefaultConeParams('cie_asano');
asanoConeParams.fieldSizeDegrees = 2*eccDeg;
asanoConeParams.ageYears = ageISETBio;
asanoConeParams.pupilDiamMM = pupilDiamMM;
asanoConeParams.indDiffParams = theData.colorDirectionParams.psiParamsStruct.coneParams.indDiffParams;
[~,~,~,~,cieConeParams,cieStaticParams] = ...
    ComputeCIEConeFundamentals(MakeItS(theData.colorDirectionParams.wls),asanoConeParams.fieldSizeDegrees,...
    asanoConeParams.ageYears,asanoConeParams.pupilDiamMM,[],[],[], ...
    [],[],[],asanoConeParams.indDiffParams);

% Photopigment absorbance
temp = load('T_log10coneabsorbance_ss');
photopigmentAbsorbance = 10.^SplineCmf(temp.S_log10coneabsorbance_ss,temp.T_log10coneabsorbance_ss,theData.colorDirectionParams.wls,2);
clear temp

% Get axial density.  This scales the normalized photopigment absorbance we
% computed just above.
axialDensity = cieConeParams.axialDensity .* (1 + asanoConeParams.indDiffParams.dphotopigment/100);

% Get OS length so we can compute axial density from first principles
%
% Rodieck, The First Steps in Seeing, p. 472 gives a value for
% axial specific density for rods and cones as about 0.015 /um.
[osLengthMicrons, osLengthMicronsFoveal] = cMosaic.outerSegmentLengthFromEccentricity(eccDeg);
specificDensity = 0.015;
axialDensityFromOSLength = specificDensity*osLengthMicrons;
fprintf('Axial density at %d degs from CIE: %0.3f, %0.3f, %0.3f; from OS length and specific density %0.3f\n', eccDeg, ...
    cieConeParams.axialDensity(1),cieConeParams.axialDensity(2),cieConeParams.axialDensity(3),axialDensityFromOSLength);

% Set up photopigment object
%
% This is based on the axialDensity (parameter name 'opticalDensity'),
% absorbance (parameter name 'absorbance'), and quantalEfficiency
% (parameter name 'peakEfficiency'). Someday we will rewrite so that naming 
% conventions are consistent across different places we implement them, perhaps.
% 
% ISETBio wants absorbance in columns, while PTB wants it in rows, so we
% also need to transposes.
photopigmentObject = cPhotoPigment('wave', theData.colorDirectionParams.wls,...
    'opticalDensity',axialDensity,'absorbance',photopigmentAbsorbance', ...
    'peakEfficiency',cieStaticParams.quantalEfficiency );

% Get the cone fundamentals we think we are using
T_cones_ISETBio = photopigmentObject.energyFundamentals';
for cc = 1:3
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) .* lensTransmittance;
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) .* macTransmittance;
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) / max(T_cones_ISETBio(cc,:));
end
fundamentalsFig = figure; clf; hold on;
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(1,:),'r','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(2,:),'g','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(3,:),'b','LineWidth',6);
xlabel('Wavelength (nm)');
ylabel('Exication Probability');
figure(fundamentalsFig);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(1,:),'k-','LineWidth',3);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(2,:),'k-','LineWidth',3);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(3,:),'k-','LineWidth',3);

%% Create and setup cone mosaic
%
% Use the ISETPipelineToolbox machinery.

%% AO version allows us to get rid of optics in effect
if (aoRender)
    if (eccVars)
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'subjectID', 0, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', 0, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end

    % Normal optics structure. Allow specified defocus.
else
    if (eccVars)
        % Build normal optics structure.
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end
end

%% Compute optical image
ss = 1; cc = 1;
theConeMosaic.PSF = oiSet(theConeMosaic.PSF,'lens',lensObject);
theOI = oiCompute(theData.ISETBioGaborObject.ISETBioGaborScene{ss,cc}, theConeMosaic.PSF);

% Some plots
sampledWavelengths = [420:40:700];
psfRangeArcMin = 10;
psfFig = visualizePSFsubStack(theOI, sampledWavelengths, psfRangeArcMin);
saveas(psfFig,fullfile(outputDir,'psfFigure.tiff'),'tiff');

uData = oiPlot(theOI, 'irradiance image with grid', [], 40);
saveas(gcf,fullfile(outputDir,'oiIrradiance.tiff'),'tiff');

%% Compute cone responses
theConeResponses = theConeMosaic.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered', ...
    'lowOpticalImageResolutionWarning',true);
%theConeMosaic.Mosaic.plot('excitations',theConeResponses);

%% Visualize the mosaic itself
mosaicFigStruct = theConeMosaic.Mosaic.visualize;
title('Cone Mosaic');
saveas(mosaicFigStruct.figureHandle,fullfile(outputDir,'mosaicImage.tiff'),'tiff');

%% Visualize noise-free responses
%
% Basic activation map
mosaicActivationFigStruct = theConeMosaic.Mosaic.visualize( ...
    'activation', theConeResponses ...
    );
title('Activation Map');
saveas(mosaicActivationFigStruct.figureHandle,fullfile(outputDir,'activationMap.tiff'),'tiff');

% Visualize a cut
horizontalSliceYcoordDegs = 0;
horizontalSliceThicknessDegs = 0.08;
coneTypesToVisualize = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
maxResponse = 5000; max(theConeResponses(:));
% visualizedResponseScalingDegs = 0.5;
% hFig = theConeMosaic.Mosaic.visualizeHorizontalConeActivationProfiles(...
%     theConeResponses, coneTypesToVisualize, ...
%     horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
%     maxResponse, visualizedResponseScalingDegs);
mosaicSliceFig = theConeMosaic.Mosaic.visualizeHorizontalConeActivationProfiles(...
    theConeResponses, coneTypesToVisualize, ...
    horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
    maxResponse,  ...
    'onTopOfConeMosaicActivation', false, ...
    'figureHandle', [], ...
    'axesHandle', [], ...
    'fontSize', 20);
saveas(mosaicSliceFig,fullfile(outputDir,'mosaicSliceFig.tiff'),'tiff');

%% Choose a circular ROI that captures most of the field and compute contrast
%
% We compute contrast across all cones of each type in the ROI.  Using the
% circular ROI eliminates edge effects.
roiCircle = regionOfInterest('shape','ellipse',...
    'center',[eccXDeg eccYDeg],...
    'majorAxisDiameter',fieldSizeDeg*0.75,...
    'minorAxisDiameter',fieldSizeDeg*0.75);

% Visualize the ROI to make sure it is in a sensible place
theConeMosaic.Mosaic.plot('roi',theConeResponses, 'roi',roiCircle);

% Get L, M, and S cone responses out of ROI.  Assumes just one
% instance/frame, which is what we are running.
lConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','L'));
mConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','M'));
sConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','S'));
lConeContrast = (max(lConeResponses(:))-min(lConeResponses(:)))/(max(lConeResponses(:))+min(lConeResponses(:)));
mConeContrast = (max(mConeResponses(:))-min(mConeResponses(:)))/(max(mConeResponses(:))+min(mConeResponses(:)));
sConeContrast = (max(sConeResponses(:))-min(sConeResponses(:)))/(max(sConeResponses(:))+min(sConeResponses(:)));

% Compute residual contrast seen by each class of cones
fprintf('L cone contrast: %0.1f%%; M cone contrast: %0.1f%%; S cone contrast: %0.1f%%\n',100*lConeContrast,100*mConeContrast,100*sConeContrast);
lmsConeContrast = [lConeContrast ; mConeContrast; sConeContrast];

%% Verify that the effective spectral sensitivities in ISETBio match intended
%
% Produce a set of monochromatic images at each sample wavelength.  The photon
% level is scene radiance in photons/sr-m2-nm-sec.
%
% This is slow but checked out well for cases I tried.
VERIFY_ISETBIO_FUNDAMENTALS = false;
if (VERIFY_ISETBIO_FUNDAMENTALS)
    pixelSize = sceneGet(theData.ISETBioGaborObject.ISETBioGaborScene{ss,cc},'rows');
    photonsPerSrM2NMSec = 1e25;
    for ww = 1:length(theData.colorDirectionParams.wls)
        % Set up a dummy scene. Spatially uniform with a black body
        % spectrum of 5000 degK.  We'll replace the scene contents just below.
        scene{ww} = sceneCreate('uniform bb',pixelSize,5000,theData.colorDirectionParams.wls);

        % Use small field of view to minimize effects of eccentricity, and also
        % so we don't need too many pixels (for efficiency in this demo).
        scene{ww} = sceneSet(scene{ww},'fov',stimulusSizeDeg);

        % Get photons and rewrite to be monochromatic constant power in
        % photons/sec-nm.
        photons = sceneGet(scene{ww},'photons');
        photons = zeros(size(photons));
        photons(:,:,ww) = photonsPerSrM2NMSec*ones(pixelSize,pixelSize);
        scene{ww} = sceneSet(scene{ww},'photons',photons);
    end

    %% Compute the retinal image and cone excitations for each scene
    %
    % Use this to get the cone fundamental that ISETBio is effectively
    % using.
    T_quantalFundamentalsFromCMosaic = zeros(3,length(theData.colorDirectionParams.wls));
    T_energyFundamentalsIFromCMosaic = zeros(3,length(theData.colorDirectionParams.wls));

    for ww = 1:length(theData.colorDirectionParams.wls)
        % Compute retinal image
        oiComputed{ww} = oiCompute(scene{ww}, theConeMosaic.PSF);

        % Compute noise free cone excitations
        coneExcitations{ww} = theConeMosaic.Mosaic.compute(oiComputed{ww},'nTrials', 1);

        % Find L, M and S cone excitations at this wavelength.  This is
        % accomplished by extracting the mean response of each cone type from
        % the excitations computed just above and diviting by the input power
        % in the scene.
        %
        % We expect the answer to be proportional to the fundamentals we
        % computed above, because we have not (yet) accounted for the geometry
        % between radiance in the scene and retinal irradiance, nor the cone
        % integration time.
        % Get L, M, and S cone responses out of ROI.  Assumes just one
        % instance/frame, which is what we are running.
        lConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','L'));
        mConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','M'));
        sConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','S'));
        T_quantalFundamentalsISETBio(1,ww) = mean(lConeResponses(:))/photonsPerSrM2NMSec;
        T_quantalFundamentalsISETBio(2,ww) = mean(mConeResponses(:))/photonsPerSrM2NMSec;
        T_quantalFundamentalsISETBio(3,ww) = mean(sConeResponses(:))/photonsPerSrM2NMSec;
    end

    for ii = 1:3
        T_energyFundamentalsIFromCMosaic(ii,:) = EnergyToQuanta(theData.colorDirectionParams.wls,T_quantalFundamentalsISETBio(ii,:)')';
        T_energyFundamentalsIFromCMosaic(ii,:) = T_energyFundamentalsIFromCMosaic(ii,:)/max(T_energyFundamentalsIFromCMosaic(ii,:));
    end

    figure(fundamentalsFig);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(1,:),'y:','LineWidth',3);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(2,:),'y:','LineWidth',3);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(3,:),'y:','LineWidth',3);
end
saveas(fundamentalsFig,fullfile(outputDir,'fundamentalsFig.tiff'),'tiff');

%% Save out what we need 
colorDirectionParams = theData.colorDirectionParams;
save(mosaicOutputFilename, ...
    'lConeResponses','mConeResponses','sConeResponses', ...
    'lConeContrast','mConeContrast','sConeContrast', 'lmsConeContrast', ...
    'T_cones_ISETBio', 'colorDirectionParams', 'eccXDeg', 'eccYDeg', 'fieldSizeDeg', ...
    'aoRender', 'noLCA', 'ageISETBio', 'dLensISETBio', 'dMacISETBio', ...
    '-v7.3');
disp('Data has been saved successfully!');
