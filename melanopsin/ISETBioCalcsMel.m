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

%% Run out a bunch of exzmples
clear;
conditionNameList = {'MelDirected1' 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2 1 2 5 10];
gaborSdDeg = 1;
stimulusSizeDeg = 4;

fieldSizeDeg = 5;
eccXDegList = [-10 15];
aoRender = false;
noLCA = false;
ageISETBio = 32;
dLensISETBio = 0;
dMacISETBio = 0;

for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        for ee = 1:length(eccXDegsList)
        ISETBioCalcsMel(conditionNameList{cc},sineFreqCyclesPerDegList(ss), ...
                gaborSdDeg,stimulusSizeDeg, ...
                fieldSizeDeg, ...
                eccXDegList(ee), ...
                aoRender, noLCA, ...
                ageISETBio, dLensISETBio, dMacISETBio ...
            );
        end
    end
end
%}

%{
%% Baseline condition test
clear;
conditionName = 'MelDirected1';
sineFreqCyclesPerDeg = 0.2;
gaborSdDeg = 1;
stimulusSizeDeg = 4;

fieldSizeDeg = 5;
eccXDeg = -10;
aoRender = false;
noLCA = false;
ageISETBio = 32;
dLensISETBio = 0;
dMacISETBio = 0;

ISETBioCalcsMel(conditionName,sineFreqCyclesPerDeg, ...
        gaborSdDeg,stimulusSizeDeg, ...
        fieldSizeDeg, ...
        eccXDeg, ...
        aoRender, noLCA, ...
        ageISETBio, dLensISETBio, dMacISETBio ...
    );
%}


function ISETBioCalcsMel(conditionName,sineFreqCyclesPerDeg, ...
            gaborSdDeg,stimulusSizeDeg, ...
            fieldSizeDeg, ...
            eccXDeg, ...
            aoRender,noLCA, ...
            ageISETBio, dLensISETBio, dMacISETBio ...
        )
% ISETBioCalcsMel

%% Initialize
close all;

%% Parameter specifications for scene that are not passed.
screenGammaMethod = 2;
sceneInputStr = sprintf('%s_Size_%0.1f_Sf_%0.1f_Sd_%0.1f_GammaMethod_%d', ...
    conditionName,stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenGammaMethod);

% Load the scene data according to parameters above
sceneInputFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
sceneInputSubdir = fullfile(sceneInputStr);
sceneInputFilename = fullfile(sceneInputFiledir,sceneInputSubdir,sprintf('sceneOutput_%s',sceneInputStr));
theData = load(sceneInputFilename);
disp('Data loaded');

%% Mosaic settings not so likely to vary
eccYDeg = 0;
subjectID = 6;
mosaicOutputStr = sprintf('%s_Size_%0.1f_Sf_%0.1f_Sd_%0.1f_GammaMethod_%d', ...
    conditionName,stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenGammaMethod);

%% Use ISETPipelineToolbox wrapper as a way to the isetbio computations
%
% These parameters define key addtional properties we are likely to vary
% Not likely to vary these in the short run, so not coded in output name for now
randSeed = false;                              % False means, no randomness
eccVars = false;
pupilDiamMM = 3;
defocusDiopters = 0;
zernikeDataBase = 'Polans2015';

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

%% Macular transmittance for the specified field size
macTransmittance = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',sqrt(eccXDeg^2 + eccYDeg^2));
unadjustedMacDensity = -log10(macTransmittance);
adjustedMacDensity = unadjustedMacDensity * (1 + dMacISETBio/100);
macTransmittance = 10.^-adjustedMacDensity;
macObject = Macular('wave',theData.colorDirectionParams.wls,'unitDensity',-log10(macTransmittance),'density',1);

%% Define a dummy mosaic that lets us get OS length at the speified eccentricity.
dummyMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'rodIntrusionAdjustedConeAperture', true, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);

%% Photopigment absorptance
%
% The individual differences model allows for shifting the absorbance along
% the wavelength axis. Do that here.
asanoConeParams = DefaultConeParams('cie_asano');
asanoConeParams.fieldSizeDegrees = sqrt(eccXDeg^2 + eccYDeg^2);
asanoConeParams.ageYears = ageISETBio;
asanoConeParams.pupilDiamMM = pupilDiamMM;
asanoConeParams.indDiffParams = theData.colorDirectionParams.psiParamsStruct.coneParams.indDiffParams;
[~,T_quantalAbsorptionProb,T_quantalExcitationProb,adjIndDiffParams,cieConeParams,cieStaticParams] = ...
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
            'rodIntrusionAdjustedConeAperture', true, ...
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
            'rodIntrusionAdjustedConeAperture', true, ...
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
            'rodIntrusionAdjustedConeAperture', true, ...
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
            'rodIntrusionAdjustedConeAperture', true, ...
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
%oiWindow(theOI);

% Compute cone responses
theConeResponses = theConeMosaic.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered', ...
    'lowOpticalImageResolutionWarning',true);
%theConeMosaic.Mosaic.plot('excitations',theConeResponses);

% Noise-free response
theConeMosaic.Mosaic.visualize( ...
    'activation', theConeResponses ...
    );

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

%% Visualize the cone response spatial profiel
% roiLine = regionOfInterest('shape', 'line', ...
%     'from', [-.25*fieldSizeDegs 0], 'to', [0.25*fieldSizeDegs,0], ...
%     'thickness', 0.1);
% theConeMosaic.Mosaic.plot('roi',theConeResponses,'roi',roiLine);

%% Save out what we need 
mosaicOutputFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
mosaicOutputFilename = fullfile(mosaicOutputFiledir,sceneInputSubdir,[sprintf('mosaicOutput_%s',mosaicOutputStr) '.mat']);
save(mosaicOutputFilename, ...
    'lConeResponses','mConeResponses','sConeResponses', ...
    'lConeContrast','mConeContrast','sConeContrast', ...
    '-v7.3');
disp('Data has been saved successfully!');
