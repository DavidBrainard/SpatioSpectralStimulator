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
%                                   gabor contrast image. It is possible to
%                                   pass multiple levels more than one.
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
%    measure                      - Default set to false. If it sets
%                                   to true, it measures the channel
%                                   primaries to calculate the point cloud.
%    verbose                      - Boolean. Default true. Controls
%                                   plotting and printout.
%    verboseDetail                - Default set to false. This prints out
%                                   all the graphs, ISETBio scene window,
%                                   and more status messages. We may not
%                                   always want to see the detailed graphs
%                                   during the process making multiple
%                                   contrast gabor images.
%
% See also: SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%           SpectralCalISETBioUsingSubroutinesV2, t_CSFGeneratorExperiment

% History:
%    01/26/22  smo        Started on it to make it as a separate function.
%    01/31/22  smo        Now you can pass multiple target gabor contrasts
%                         to generate the images at once.

%% Set parameters.
arguments
    targetContrast {mustBeInRange(targetContrast,0,1,'inclusive')}
    colorDirectionParams
    spatialTemporalParams
    options.measure (1,1) = false
    options.verbose (1,1) = true
    options.verboseDetail (1,1) = false
end

%% Say hello.
if (options.verbose)
    fprintf('Starting to make gabor image(s)...\n');
end

%% Set the target contrast gabor here.
%
% This part is for putting multiple contrast levels in the params
% structure.
colorDirectionParams.spatialGaborTargetContrast = targetContrast;

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

%% Do all calibraiton loading.
screenGammaMethod = 2;
[screenCalObj,channelCalObjs] = LoadAndSetExperimentCalFiles(colorDirectionParams,'screenGammaMethod',screenGammaMethod,'verbose',options.verboseDetail);

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
[screenPrimaryChannelObject,backgroundChannelObject] = SetupChannelPrimaries(colorDirectionParams,channelCalObjs,projectIndices,'verbose',options.verboseDetail);

%% Set the screen primaries.
%
% We want these to match those we set up with the channel calculations
% above.  Need to reset sensor color space after we do this, so that the
% conversion matrix is properly recomputed.
screenCalObj.set('P_device', screenPrimaryChannelObject.screenPrimarySpd);
SetSensorColorSpace(screenCalObj, colorDirectionParams.T_cones, colorDirectionParams.S);

%% Create ISETBio display from the calibration file.
[ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj,'verbose',options.verboseDetail);

%% Set up the background screen primaries.
backgroundScreenPrimaryObject = SetupBackground(colorDirectionParams,screenCalObj,backgroundChannelObject,'verbose',options.verboseDetail);

%% Make a monochrome Gabor patch in range -1 to 1.
%
% This is our monochrome contrast modulation image. Multiply by the max
% contrast vector to get the LMS contrast image. The function includes the
% quantization of the gabor image.
nQuantizeBits = 14;
[rawMonochromeUnquantizedContrastGaborImage, rawMonochromeUnquantizedContrastGaborCal, rawMonochromeContrastGaborCal, ...
    stimulusN, centerN, stimulusHorizSizeDeg, stimulusHorizSizeMeters] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,'verbose',options.verboseDetail,'nQuantizeBits',nQuantizeBits);

%% Get cone contrast/excitation gabor image.
[ptCldObject,standardGaborCalObject,screenCalObj,backgroundScreenPrimaryObject] = ...
    SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,...
    screenCalObj,backgroundScreenPrimaryObject,screenPrimaryChannelObject,...
    'measure',options.measure,'warmupTimeMinutes',30,'verbose',options.verboseDetail);

%% Make image from point cloud.
gaborImageObject = MakeImageSettingsFromPtCld(ptCldObject,screenCalObj,standardGaborCalObject,...
    backgroundScreenPrimaryObject.screenBgExcitations,stimulusN,'verbose',options.verboseDetail);

%% Put the image into an ISETBio scene.
ISETBioGaborObject = MakeISETBioSceneFromImage(colorDirectionParams,gaborImageObject,standardGaborCalObject,...
    ISETBioDisplayObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg,'verbose',options.verboseDetail);

%% Save out the images in a single variable.
gaborRGBImage = gaborImageObject.standardSettingsGaborImage;
gaborISETBioScene = ISETBioGaborObject.ISETBioGaborScene;

%% Say goodbye. 
if (options.verbose)
    fprintf('Gabor image(s) has(have) been successfully created! \n');
end

end
