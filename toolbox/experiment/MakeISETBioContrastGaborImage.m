function [gaborISETBioScene,gaborRGBImage] = MakeISETBioContrastGaborImage(targetContrast,colorDirectionParams,spatialTemporalParams,options)
% Make a contrast gabor image in both image and ISETBio scene formats.
%
% Syntax:
%    [gaborISETBioScene,gaborRGBImage] = MakeISETBioContrastGaborImage(targetContrast,spatialTemporalParams)
%
% Description:
%    This makes a contrast gabor image in both image and ISETBio scene
%    formats. This used to be a long script (SpectralCalISETBio) and this
%    is the version of a single function to create an image with different
%    contrast levels.
%
%    This is basically the function version of the
%    (SpectralCalISETBioUsingSubroutinesV2).
% 
% Inputs:
%    targetContrast               - Desired maximum contrast to have in the
%                                   gabor contrast image.
%    colorDirectionParams         - Structure with the parameters to
%                                   calculate a contrast gabor image.
%    spatialTemporalParams        - Structure with the parameters to decide
%                                   the spatial temporal characteristics of
%                                   the gabor image.
%
% Outputs:
%    gaborISETBioScene            - Created gabor image in a ISETBio scene
%                                   format.
%    gaborRGBImage                - Created gabor image in RGB format.
%
% Optional key/value pairs:
%    verbose                      - Boolean. Default true. Controls
%                                   plotting and printout.
%
% See also: SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%           SpectralCalISETBioUsingSubroutinesV2, t_CSFGeneratorExperiment

% History:
%    01/26/22  smo        Started on it to make it as a separate function.

%% Set parameters.
arguments
    targetContrast (1,1) {mustBeInRange(targetContrast,0,1,'inclusive')}
    colorDirectionParams
    spatialTemporalParams
    options.verbose (1,1) = true
end

%% Say hello.
if (options.verbose)
    fprintf('Starting to create a Gabor image with the contast (%2.2f)...\n',targetContrast);
end

%% Check if we have all spatial temporal parameters that we need.
spatialTemporalParamsCheck = {'sineFreqCyclesPerDeg','gaborSdDeg','stimulusSizeDeg'};
nspatialTemporalParamsCheck = size(spatialTemporalParamsCheck,2);
for cc = 1:nspatialTemporalParamsCheck
    if(any(~isfield(spatialTemporalParams,spatialTemporalParamsCheck{cc})))
        error(append('Not enough spatial temporal parameters, missing parameter: ',spatialTemporalParamsCheck{cc}));
    end
end

% Set image spatial parameters here.
sineFreqCyclesPerDeg = spatialTemporalParams.sineFreqCyclesPerDeg;
gaborSdDeg = spatialTemporalParams.gaborSdDeg;
stimulusSizeDeg = spatialTemporalParams.stimulusSizeDeg;

%% Verbose setting for printing out all graphs.
%
% Set to true to get more output. This is different one than
% options.verbose. If it sets to true, it prints out all the graphs
% including the scene image on ISETBio, but we don't want to print out all
% that stuffs when we run the experiment.
VERBOSEDETAIL = false;

%% Do all calibraiton loading.
screenGammaMethod = 2;
[screenCalObj,channelCalObjs] = LoadAndSetExperimentCalFiles(colorDirectionParams,'screenGammaMethod',screenGammaMethod,'verbose',VERBOSEDETAIL);

%% Use extant machinery to get primaries from spectrum.
%
% Define wavelength range that will be used to enforce the smoothness
% through the projection onto an underlying basis set.  We don't the whole
% visible spectrum as putting weights on the extrema where people are not
% sensitive costs us smoothness in the spectral region we care most about.
lowProjectWl = 400;
highProjectWl = 700;
projectIndices = find(colorDirectionParams.wls > lowProjectWl & colorDirectionParams.wls < highProjectWl);

%% Find primaries with desired LMS contrast.
[screenPrimaryChannelObject,backgroundChannelObject] = SetupChannelPrimaries(colorDirectionParams,channelCalObjs,projectIndices,'verbose',VERBOSEDETAIL);

%% Set the screen primaries.
%
% We want these to match those we set up with the channel calculations
% above.  Need to reset sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device', screenPrimaryChannelObject.screenPrimarySpd);
SetSensorColorSpace(screenCalObj, colorDirectionParams.T_cones, colorDirectionParams.S);

%% Create ISETBio display from the calibration file.
[ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj,'verbose',VERBOSEDETAIL);

%% Set up the background screen primaries.
backgroundScreenPrimaryObject = SetupBackground(colorDirectionParams,screenCalObj,backgroundChannelObject,'verbose',VERBOSEDETAIL);

%% Make a monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image. Multiply by the max
% contrast vector to get the LMS contrast image. The function includes the
% quantization of the gabor image.
nQuantizeBits = 14;
[rawMonochromeUnquantizedContrastGaborImage, rawMonochromeUnquantizedContrastGaborCal, rawMonochromeContrastGaborCal, ...
    stimulusN, centerN, stimulusHorizSizeDeg, stimulusHorizSizeMeters] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,'verbose',VERBOSEDETAIL,'nQuantizeBits',nQuantizeBits);

%% Get cone contrast/excitation gabor image.
[ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,...
    screenCalObj,backgroundScreenPrimaryObject.screenBgExcitations,'verbose',VERBOSEDETAIL);

%% Make image from point cloud.
gaborImageObject = MakeImageSettingsFromPtCld(ptCldObject,screenCalObj,standardGaborCalObject,...
    backgroundScreenPrimaryObject.screenBgExcitations,stimulusN,'verbose',VERBOSEDETAIL);

%% Put the image into an ISETBio scene.
ISETBioGaborObject = MakeISETBioSceneFromImage(colorDirectionParams,gaborImageObject,standardGaborCalObject,...
    ISETBioDisplayObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg,'verbose',VERBOSEDETAIL);

%% Save out the images in a single variable.
gaborRGBImage = gaborImageObject.standardSettingsGaborImage;
gaborISETBioScene = ISETBioGaborObject.ISETBioGaborScene;

%% Say goodbye. 
%
% Print out if everything goes well.
if (options.verbose)
    fprintf('Gabor image with the contast (%2.2f) has been successfully created! \n',targetContrast);
end

end
