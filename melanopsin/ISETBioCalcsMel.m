% ISETBioCalcsMel

% Initialize
clear; close all;

% Load the data
testFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
testFilename = fullfile(testFiledir,'testImageDataISETBio');
theData = load(testFilename);
disp('Data loaded');

% Use recon pipeline wrapper as a way to the isetbio computations, because
% it is nicely encapsulated.
%
% These parameters define the mosaic
randSeed = false;
pupilDiamMM = 3;
aoRender = false;
noLCA = false;
defocusDiopters = 0;
eccVars = false;
fieldSizeDegs = 5;
circleSizeDegs = 3;
eccXDegs = 9;
eccYDegs = 0;
addPoissonNoise = true;
subjectID = 6;
zernikeDataBase = 'Polans2015';

if (eccVars)
    % Build normal optics structure.
    theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
        'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
        'subjectID', subjectID, ...
        'noLCA', noLCA, ...
        'zernikeDataBase', zernikeDataBase);
else
    theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
        'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
        'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
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

% Compute optical image
ss = 1; cc = 1;
theOI = oiCompute(theData.ISETBioGaborObject.ISETBioGaborScene{ss,cc}, theConeMosaic.PSF);
oiWindow(theOI);


% Compute cone responses
theConeResponses = theConeMosaic.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered', ...
    'lowOpticalImageResolutionWarning',true);
%theConeMosaic.Mosaic.plot('excitations',theConeResponses);

% % Compute residual contrast seen by each class of cones
% lConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.lConeIndices));
% lConeContrast = (max(lConeResponses(:))-min(lConeResponses(:)))/(max(lConeResponses(:))+min(lConeResponses(:)));
% mConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.mConeIndices));
% mConeContrast = (max(mConeResponses(:))-min(mConeResponses(:)))/(max(mConeResponses(:))+min(mConeResponses(:)));
% sConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.sConeIndices));
% sConeContrast = (max(sConeResponses(:))-min(sConeResponses(:)))/(max(sConeResponses(:))+min(sConeResponses(:)));
% fprintf('L cone contrast: %0.1f%%; M cone contrast: %0.1f%%; S cone contrast: %0.1f%%\n',100*lConeContrast,100*mConeContrast,100*sConeContrast);

% Noise-free response
theConeMosaic.Mosaic.visualize( ...
    'activation', theConeResponses ...
    );

% Now choose a small circular ROI
roiCircle = regionOfInterest('shape','ellipse',...
    'center',[eccXDegs eccYDegs],...
    'majorAxisDiameter',fieldSizeDegs*0.5,...
    'minorAxisDiameter',fieldSizeDegs*0.5);

theConeMosaic.Mosaic.plot('roi',theConeResponses, 'roi',roiCircle);

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