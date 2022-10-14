% t_CSFGeneratorMakeTestImages.
%
% This makes contrast gabor images for the SACC experiment.
%
% See Also:
%    t_CSFGeneratorExperiment.

% History:
%    10/12/22   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set paramters.
sineFreqCyclesPerDegTarget = [3 6 9 12 18];
minContrast = 0.0003;
nContrasts = 30;
measure = true;
conditionName = 'LminusMSmooth';

%% Make a loop here to make all test images at once for all spatial frequencies.
for ss = 1:numel(sineFreqCyclesPerDegTarget)
    sineFreqCyclesPerDeg = sineFreqCyclesPerDegTarget(ss);
    
    % Set up color direction
    %
    % Set spatialGaborTargetContrast differently according to spatial frequency.
    if (sineFreqCyclesPerDeg == 3)
        spatialGaborTargetContrast = 0.02;
    elseif (sineFreqCyclesPerDeg == 18)
        spatialGaborTargetContrast = 0.06;
    else % This is for 6, 9, and 12 cpd.
        spatialGaborTargetContrast = 0.04;
    end
    
    colorDirectionParams = SetupColorDirection(conditionName,...
        'spatialGaborTargetContrast',spatialGaborTargetContrast);
    
    %% Image spatial parameters.
    %
    % Control sineFreqCyclesPerDeg for changing spatial frequncy of the
    % contrast gabor pattern. For example, setting it to 1 means 1 cpd.
    %
    % We used to make the contrast gabor size in 1.5 gaborSdDeg, and now we
    % are trying to test the size of 0.75.
    spatialTemporalParams.sineFreqCyclesPerDeg = sineFreqCyclesPerDeg;
    spatialTemporalParams.gaborSdDeg = 0.75;
    spatialTemporalParams.stimulusSizeDeg = 7;
    spatialTemporalParams.sineImagePhaseShiftDeg = [0];
    
    %% Instantiate a sceneEngine.
    %
    % First set up the scene parameters that will be needed by
    % the sceSACCDisplay.
    %
    % First step is to predefine the contrasts that we will allow the
    % psychophysics to work over. This gives us a finite list of scenes
    % to compute for.
    experimentParams.minContrast = minContrast;
    experimentParams.nContrasts = nContrasts;
    experimentParams.stimContrastsToTest = [0 logspace(log10(experimentParams.minContrast), ...
        log10(colorDirectionParams.spatialGaborTargetContrast), experimentParams.nContrasts)];
    
    experimentParams.slopeRangeLow = 0.5;
    experimentParams.slopeRangeHigh = 6;
    experimentParams.slopeDelta = 0.5;
    experimentParams.nTest = 1;
    experimentParams.nQUESTEstimator = 1;
    
    % If this is set to true, we will measure primaries and use them.
    experimentParams.measure = measure;
    
    %% Make contrast gabor images and save.
    % Now do all the computation to get us ISETBio scenes and RGB images for
    % each predefined contrast, relative to the parameters set up above.
    %
    % Move the precomputed data into the format for the sceSACCDisplay scene
    % engine. This will take some fair amount of time to run it.
    %
    % For experiment, turn off generation of ISETBio scenes because it eats up
    % time and even worse memory.  But can turn on in the future for
    % computational analyses.
    %
    % Also, lightVer prints out less variables inside the function so that we
    % can save a lot of memory and time. What we get is the same.
    noISETBio = true;
    lightVer = true;
    printGaborSpds = false;
    
    % Make contrast gabor images here.
    [sceneParamsStruct.predefinedSceneSequences, sceneParamsStruct.predefinedRGBImages, experimentParams.screenPrimarySettings desiredSpdGaborCal] = ...
        MakeISETBioContrastGaborImage(experimentParams.stimContrastsToTest, ...
        colorDirectionParams,spatialTemporalParams,'measure',experimentParams.measure,...
        'verbose',true,'noISETBio',noISETBio,'lightVer',lightVer,'printGaborSpds',printGaborSpds);
    
    % Set some of the scene parameters.
    sceneParamsStruct.predefinedContrasts = experimentParams.stimContrastsToTest;
    
    % Save the images and params.
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,'TestImages',sprintf('RunExpData_%d_cpd_%s.mat',...
            sineFreqCyclesPerDeg,dayTimestr));
        save(testFilename,'colorDirectionParams','spatialTemporalParams','sceneParamsStruct', ...
            'experimentParams','noISETBio','lightVer');
    end
    
    % Show the progress.
    fprintf('\t Progress making gabor image - Spatial frequency (%d/%d) \n', ss, numel(sineFreqCyclesPerDegTarget));
end

%% Move the measured screen primaries data in the different folder.
%
% Here we will save the measurement data in the separate folder and delete
% the temporary saved data.
STOREMEASUREMENT = true;

if (STOREMEASUREMENT)
    if (ispref('SpatioSpectralStimulator','SACCData'))
        % Load the measurement file name.
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        testFilename = fullfile(testFiledir,'TestImages','MeasurementData','targetScreenSpdMeasured.mat');
        load(testFilename);
        
        % Save it with the date.
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFiledirNew = fullfile(testFiledir,'TestImages','MeasurementData');
        testFilenameNew = fullfile(testFiledirNew,sprintf('targetScreenSpdMeasured_%s.mat',dayTimestr));
        save(testFilenameNew,'targetScreenSpdMeasured');
        
        % Delete the old file.
        delete(testFilename);
    end
end