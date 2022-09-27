% Show how to use CSF generator engine to run an experiment.
%
% Syntax:
%    t_CSFGeneratorExperiment
%
% Description:
%    Demonstrates how to use our CSF generator machinery to drive a
%    psychophysical experiment.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    ISETBioCSFGenerator toolbox and its tutorials.
%
% History:
%    01/18/22  dhb, smo  - Start writing.
%    01/26/22  smo       - Added the part making the contrast gabor
%                          image with different contrast levels.
%    02/02/22  smo       - Updated on the function
%                          computePerformanceSACCDisplay for getting
%                          response from patients.
%    02/08/22  dhb,smo   - Added an option to skip making the ISETBio scenes
%                          which takes time and memory a lot.
%    03/17/22  smo       - Now we set the temporal scene support time
%                          differently for stimuli and cross fixation point.
%    04/27/22  smo       - Added a PTB running option for displaying test
%                          gabor image in either vertical or horizontal
%                          direction.
%    05/09/22  smo       - Added an option to make phase shift on gabor
%                          image.
%    06/23/22  smo       - Now we save the gabor image data (RunExpData) by
%                          containing its spatial frequency info in its
%                          file name.
%    07/11/22  smo       - Added an option to display stimuli in debug mode
%                          to keep displaying images until button pressed.
%    07/12/22  smo       - Now we can print out the spds of gabor images
%                          which is useful to create sRGB images if we
%                          want.
%    07/13/22  smo       - Added an option to make stimuli
%                          presentation gradually ramping on and off.
%    07/18/22  smo       - Added an option to make a time delay on null
%                          image before showing the test contrast image.
%    07/20/22  smo       - Added a practice trials before the main
%                          experiment.
%    08/08/22  smo       - Added an option to display a small focusing
%                          point on the test contrast images to minimize
%                          the artifacts.
%    08/29/22  smo       - Added an option to rotate the images.

%% Initialization.
clear; close all;

%% Set up parameters for making gabor images.
%
% You can load the data if you saved the images. We will load the images when
% we run the main experiment to save the time for making the images.
LOADDATA = true;
METHODOFADJUSTMENT = true;
PRACTICETRIALS = true;
conditionName = 'LminusMSmooth';

if (LOADDATA)
    % Set the condition of the images.
    sineFreqCyclesPerDeg = 3;
    gaborSdDeg = 0.75;
    SAVETHERESULTS = true;
    
    % Load the data here.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = fullfile(testFiledir,sprintf('RunExpData_%s_%d_cpd_%.2f_SdDeg.mat',...
            conditionName,sineFreqCyclesPerDeg,gaborSdDeg));
        load(testFilename);
    end
    
elseif (~LOADDATA)
    %% Set up color direction
    %
    % Set spatialGaborTargetContrast = 0.03 for the spatial frequency 18
    % cpd. If set 0.02, it's almost impossible to detect the stimuli.
    %
    % Now we set it as double.
    % 1 cpd = 0.02
    % the others = 0.04
    % 18 cpd = 0.06
    spatialGaborTargetContrast = 0.06;
    colorDirectionParams = SetupColorDirection(conditionName,...
        'spatialGaborTargetContrast',spatialGaborTargetContrast);
    
    %% Image spatial parameters.
    %
    % Control sineFreqCyclesPerDeg for changing spatial frequncy of the
    % contrast gabor pattern. For example, setting it to 1 means 1 cpd.
    %
    % We used to make the contrast gabor size in 1.5 gaborSdDeg, and now we
    % are trying to test the size of 0.75.
    spatialTemporalParams.sineFreqCyclesPerDeg = 18;
    spatialTemporalParams.gaborSdDeg = 0.75;
    spatialTemporalParams.stimulusSizeDeg = 7;
    spatialTemporalParams.sineImagePhaseShiftDeg = [0 90 180 270];
        
    %% Instantiate a sceneEngine.
    %
    % First set up the scene parameters that will be needed by
    % the sceSACCDisplay.
    %
    % First step is to predefine the contrasts that we will allow the
    % psychophysics to work over. This gives us a finite list of scenes
    % to compute for.
    experimentParams.minContrast = 0.0003;
    experimentParams.nContrasts = 20;
    experimentParams.stimContrastsToTest = [0 logspace(log10(experimentParams.minContrast), ...
        log10(colorDirectionParams.spatialGaborTargetContrast), experimentParams.nContrasts)];
    
    experimentParams.slopeRangeLow = 0.5;
    experimentParams.slopeRangeHigh = 6;
    experimentParams.slopeDelta = 0.5;
    experimentParams.nTest = 1;
    experimentParams.nQUESTEstimator = 1;
    
    % If this is set to true, we will measure primaries and use them.
    experimentParams.measure = true;
end

%% Set up experimental parameters here.
%
% Set the number of trials here.
%
% RunningMode can be chosen among three
% [PTB-sequential; PTB-directional; simulation].
experimentParams.minTrial = 40;
experimentParams.maxTrial = 40;
experimentParams.nTestValidation = 20;
experimentParams.runningMode = 'PTB-directional';
experimentParams.expKeyType = 'gamepad';
experimentParams.beepSound = true;
experimentParams.autoResponse = false;
experimentParams.debugMode = false;
experimentParams.preStimuliDelaySec = 0;
experimentParams.movieStimuli = true;
experimentParams.movieImageDelaySec = 0.2;

% Set the presentation time for each target and crossbar images in seconds
% unit.
sceneParamsStruct.predefinedTemporalSupport = 0.4;
sceneParamsStruct.predefinedTemporalSupportCrossbar = 1.0;
sceneParamsStruct.sineImagePhaseShiftDeg = spatialTemporalParams.sineImagePhaseShiftDeg;
sceneParamsStruct.addFixationPointImage = true;
sceneParamsStruct.addNoiseToImage = true;
sceneParamsStruct.rotateImageDeg = 45;

% Set numbers when using auto response.
autoResponseParams.psiFunc = @qpPFWeibullLog;
autoResponseParams.thresh = 0.004;
autoResponseParams.slope = 2;
autoResponseParams.guess = 0.5;
autoResponseParams.lapse = 0.01;
autoResponseParams.psiParams = [log10(autoResponseParams.thresh) autoResponseParams.slope autoResponseParams.guess autoResponseParams.lapse];

%% Make contrast gabor images and save.
if (~LOADDATA)
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
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = fullfile(testFiledir,sprintf('RunExpData_%s_%d_cpd_%.2f_SdDeg.mat',...
            conditionName,spatialTemporalParams.sineFreqCyclesPerDeg,spatialTemporalParams.gaborSdDeg));
        save(testFilename,'colorDirectionParams','spatialTemporalParams','sceneParamsStruct', ...
            'experimentParams','noISETBio','lightVer');
    end
end

%% Match the contrast range scale.
%
% Round the contrast to match the scale.
experimentParams.stimContrastsToTest = round(experimentParams.stimContrastsToTest,4);

% Update the scene parameters.
sceneParamsStruct.predefinedContrasts = round(experimentParams.stimContrastsToTest,4);
  
%% Open projector and set the screen primary settings as we found.
if (or(strcmp(experimentParams.runningMode,'PTB-sequential'),strcmp(experimentParams.runningMode,'PTB-directional')))
    initialScreenSetting = [0 0 0]';
    [window windowRect] = OpenPlainScreen(initialScreenSetting);
    SetChannelSettings(experimentParams.screenPrimarySettings);
    
elseif (strcmp(experimentParams.runningMode,'simulation'))
    % Set arbitrary numbers to pass and these will not be used in the
    % further function.
    window = 1;
    windowRect = [0 0 1920 1080];
end

%% Method of adjustment.
if (METHODOFADJUSTMENT)
    % Method of adjustment happens here and get the contrast range.
    [estDomainValidation estDomainValidationNominalLinear contrastFoundLinear thresholdEstLinear] = ...
        GetContrastRangeTrials(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect);
end 

% Set the auto response params empty if it is not used in the main
% experiment.
if (~experimentParams.autoResponse)
    autoResponseParams = [];
end

%% Create the scene engine.
theSceneEngine = sceneEngine(@sceSACCDisplay,sceneParamsStruct);

%% Construct a QUEST threshold estimator estimate threshold on log contrast
%
% We have found that computing out psychometric functions takes a long time
% if you use the method of constant stimuli.  Thus we use an adaptive
% psychophysical procedure, QUEST+, to make the computations efficient.
% You can learn more about QUEST+ at the gitHub site that has our QUEST+
% Matlab implementation: https://github.com/BrainardLab/mQUESTPlus.
%
% A features of QUEST+ is that you need to specify a discrete list of
% contrast values that will be tested, which in the wrapper here is also
% taken as the set of possible thresholds as QUEST+ chooses stimulus values.
% We also discretize the possible slopes of the underlying Weibull
% psychometric function.  As noted above, appropriate choices for these
% depend on the neural engine (and possibly the quality of the response
% classifier) and need to be set on the basis of experience.  There is also
% a little art to choosing the spacing of the discrete values.  Too many
% and QUEST+ will run slowly; too few and it may not be able to put stimuli
% in the right places.
estDomain  = log10(experimentParams.stimContrastsToTest(2:end));
slopeRange = experimentParams.slopeRangeLow: experimentParams.slopeDelta : experimentParams.slopeRangeHigh;

% Note the explicit setting of the PF for the questThresholdEngine.  Using
% @qpPFWeibullLog causes it all to happen in log10 units, rather than the
% default dB units.
%
% Run 'numEstimator > 1' interleaved QUEST+ objects. In this case,
% the threshold engine calculates the running standard error (SE)
% among those objects. Essentially, this is an internal estimate of
% the precision of the c10urrent threshold estimate.  Given this, the
% stopping criterion is triggered when either of the following is
% true
%   1) total number of trials >= 'minTrial' AND the 'stopCrition'
%   is true.
%   2) when total number of trials >= 'maxTrial'.
% The total number of trials includes all the trials across the
% multiple QUEST+ objects.
%
% The 'stopCriterion' itself can be defined in two ways.
%   1) A single number. In this case, the criterion will be
%      true when the SE is less than the specified number.
%   2) A function handle that takes the current estimate of threshold
%      and SE estimates as input arguments, and returns a boolean variable.
%      For example: we can use a relative criterion w.r.t. the magnitude of threshold
%
% Comment in either of the lines below to chose which you do. As
% with the contrast range, there is some art to choosing the
% specific values for the stop criterion, as well as the min and
% max number of trials.  For example, the relative stopping
% criterion in the function handle below can terminate too early if
% initial threshold values are large.  This can be avoided by
% appropriate choice of minimum number of trials.
experimentMode = 'validation';

switch experimentMode
    case 'adaptive'
        stopCriterion = @(threshold, se) se / abs(threshold) < 0.01;
        % Set up the estimator object.
        estimator = questThresholdEngine('minTrial', experimentParams.minTrial, ...
            'maxTrial', experimentParams.maxTrial, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, ...
            'numEstimator', experimentParams.nQUESTEstimator, ...
            'stopCriterion', stopCriterion, 'qpPF',@qpPFWeibullLog);
        
    case 'validation'
        % Set the test contrast domain to validate.
        %
        % We set the contrast range differently according to spatial
        % frequencies. These ranges are with the size of pupil 3.0 mm,
        % which works fine for Semin and David.
        if (~METHODOFADJUSTMENT)
            switch sineFreqCyclesPerDeg
                case 3
                    lowerLimEstDomain  = 0.0019;
                    higherLimEstDomain = 0.0046;
                    
                    % Reset it for Geoff.
                    lowerLimEstDomain  = 0.0035;
                    higherLimEstDomain = 0.0062;
                case 6
                    lowerLimEstDomain  = 0.0027;
                    higherLimEstDomain = 0.0081;
                case 9
                    lowerLimEstDomain  = 0.0038;
                    higherLimEstDomain = 0.0092;
                case 12
                    lowerLimEstDomain  = 0.0038;
                    higherLimEstDomain = 0.0092;
                case 18
                    lowerLimEstDomain  = 0.0071;
                    higherLimEstDomain = 0.0152;
                    
                    % Reset it for David
                    lowerLimEstDomain  = 0.0038;
                    higherLimEstDomain = 0.0120;
                otherwise
            end
            
            % Set the contrast range here.
            estDomainIndex = find(and(experimentParams.stimContrastsToTest >= lowerLimEstDomain, ...
                experimentParams.stimContrastsToTest <= higherLimEstDomain));
            estDomainValidation = estDomain(estDomainIndex-1);
        end
        
        % Estimation domain in linear unit. 
        estDomainValidationLinear = 10.^estDomainValidation;
        
        % Set up the estimator object.
        estimator = questThresholdEngine('validation',true, ...
            'nRepeat',experimentParams.nTestValidation, 'estDomain', estDomainValidation, ...
            'slopeRange', slopeRange, 'qpPF', @qpPFWeibullLog);
        
    otherwise
        error('Experiment mode should be set either (adaptive) or (validation).');
end

%% Generate the NULL scene sequence
%
% Threshold will be measured as a perturbation from this scene.  Typically
% it will correspond to zero contrast, but it doesn't have to.
%
% This call illustrates the compute method of the sceneEngine class.  We
% pass the desired contrast, and back comes a cell array as the first
% returned value, with one ISETBio scene for each time point.  The
% corresponding times are in the second returned value, an array.  Times
% are in seconds.
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds, nullStatusReportStruct] ...
    = theSceneEngine.compute(nullContrast);
theCrossbarTemporalSupportSeconds = sceneParamsStruct.predefinedTemporalSupportCrossbar;

if (noISETBio)
    nullStatusReportStruct.RGBimage = sceneParamsStruct.predefinedRGBImages{1,1};
end

%% Threshold estimation with QUEST+.
%
% There is some logic here to cache scenes and classifiers for each
% contrast tested, so that things will run faster if a contrast is repeated
% a second time.  This can use up space, so for real problems you may need
% to think about space-time tradeoffs, caching things to disk, etc.
%
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

%% Practice trials before the main experiment if you want.
%
% Set the images to use for practice trials.
if (PRACTICETRIALS)
    
    %% Make images of the practice trials.
    %
    % We make initial, finishing images for disply and load the contrast
    % image for practice trials.
    %
    % Note that our DMD displays the vertical symmetry image of the
    % original images, so here we make vertical symmetry image to look fine
    % on the DMD.
    %
    % Make initial screen image.
    imageSize = size(nullStatusReportStruct.RGBimage,2);
    messageInitialRGBImage_1stLine = 'Press any button to start';
    messageInitialRGBImage_2ndLine = 'Main Experiment';
    initialRGBImagePractice = insertText(nullStatusReportStruct.RGBimage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialRGBImage_1stLine messageInitialRGBImage_2ndLine},...
        'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
    initialRGBImagePractice = fliplr(initialRGBImagePractice);
    
    % Set the contrast image for practice.
    practiceTestContrast = max(sceneParamsStruct.predefinedContrasts);
    practiceTestContrastIndex = find(sceneParamsStruct.predefinedContrasts == practiceTestContrast);
    practiceRGBImage = sceneParamsStruct.predefinedRGBImages{1,practiceTestContrastIndex};
    
    %% Display the initial screen of the practice trials.
    SetScreenImage(initialRGBImagePractice, window, windowRect,'verbose',true);
    
    % Press any button to proceed.
    GetGamepadResp;
    disp('Practice trial is going to be started!');
    
    %% Display the null image.
    addFixationPointCircle = 'circle';
    [imageTextureNull imageWindowRect] = MakeImageTexture(nullStatusReportStruct.RGBimage, window, windowRect, ...
    'addNoiseToImage', sceneParamsStruct.addNoiseToImage, 'addFixationPoint', addFixationPointCircle, 'verbose', false);
    FlipImageTexture(imageTextureNull, window, imageWindowRect,'verbose', false);
    
    % Press any button to proceed.
    GetGamepadResp;
    disp('Experiment is going to be started!');
    
    %% Start the practice trials here.
    %
    % Here we will display the highest contrast image five times. It
    % should be easily visible. This is for subjects to get used to how
    % they evaluate stimuli and to expect what's going to happen during the
    % experiment. Here we set 5 times, but we can increase the number of
    % the trials if needed.
    nPracticeTrials = 5;
    for pp = 1:nPracticeTrials
        
        % Print out the progress.
        fprintf('Starting practice trial (%d/%d) \n', pp, nPracticeTrials);
        
        [correct] = computePerformanceSACCDisplay(nullStatusReportStruct.RGBimage, practiceRGBImage, ...
            theSceneTemporalSupportSeconds,theCrossbarTemporalSupportSeconds,practiceTestContrast,window,windowRect,...
            'runningMode',experimentParams.runningMode,'autoResponse',autoResponseParams,...
            'expKeyType',experimentParams.expKeyType,'beepSound',experimentParams.beepSound,...
            'debugMode',experimentParams.debugMode,'movieStimuli',experimentParams.movieStimuli,...
            'movieImageDelaySec',experimentParams.movieImageDelaySec,...
            'preStimuliDelaySec',experimentParams.preStimuliDelaySec, 'addNoiseToImage', sceneParamsStruct.addNoiseToImage, ...
            'addFixationPointImage', sceneParamsStruct.addFixationPointImage,...
            'rotateImageDeg',sceneParamsStruct.rotateImageDeg, 'verbose',true);
    end
    
    disp('Practice trial has been ended!');
end

%% Main experiment starts from here.
%
% If PTB mode and not simulating response, wait for subject
% to press a button before starting trials. In the end,
% probably want to move this to a place where the test cross
% is displayed for the first time.  Or display test cross here,
% or something.
if (~PRACTICETRIALS)
    
    % Display crossbar image.
    [imageTextureNull, imageWindowRectNull] = MakeImageTexture(nullStatusReportStruct.RGBimage, window, windowRect,...
        'addFixationPoint', 'circle', 'verbose', true);
    FlipImageTexture(imageTextureNull, window, imageWindowRectNull);
    
    % Press any button to proceed.
    GetGamepadResp;
end

% All experimental trials happen here that include displaying test images
% and getting responses.
%
% And we are collecting flip time whenever the displaying image changes.
flipTime = [];
rngVal = {};

while (nextFlag)
    % Convert log contrast -> contrast.
    %
    % We round it to set it as exact the number of the target contrast.
    % Otherwise, it throws the error.
    testContrast = round(10 ^ logContrast,4);
    if (isempty(find(testContrast == experimentParams.stimContrastsToTest)))
        error('Test contrast not in predefined list. Check numerical precision');
    end
    
    % Get the scene sequence and RGB info for the desired contrast.
    % Our scene engine provides us with the RGB values we need in the
    % status report structure.  A bit of overloading of the original
    % intent, but should work just fine.
    [theTestSceneSequences, ~, testStatusReportStruct] = ...
        theSceneEngine.compute(testContrast);
    
    % Run the trial and get the response. This routine
    % takes the RGB image info for the trial and returns 1 if is a
    % correct trial and 0 if incorrect trial.
    %
    % Current version of 'computePerformanceSACCDisplay' does not use
    % displayControlStruct, which was needed in the previous version. Maybe
    % we can bring it back when it is needed.
    %
    % Get a response here. Make a loop for the number of trials.
    for tt = 1:experimentParams.nTest
        [correct(tt) flipTimeTemp(:,tt) rngValTemp{tt}] = computePerformanceSACCDisplay(...
            nullStatusReportStruct.RGBimage, testStatusReportStruct.RGBimage, ...
            theSceneTemporalSupportSeconds,theCrossbarTemporalSupportSeconds,testContrast,window,windowRect,...
            'runningMode',experimentParams.runningMode,'autoResponse',autoResponseParams,...
            'expKeyType',experimentParams.expKeyType,'beepSound',experimentParams.beepSound,...
            'debugMode',experimentParams.debugMode,'movieStimuli',experimentParams.movieStimuli,...
            'movieImageDelaySec',experimentParams.movieImageDelaySec,...
            'preStimuliDelaySec',experimentParams.preStimuliDelaySec,'addNoiseToImage',sceneParamsStruct.addNoiseToImage, ...
            'addFixationPointImage', sceneParamsStruct.addFixationPointImage, 'rotateImageDeg', sceneParamsStruct.rotateImageDeg,'verbose',true);
    end
    
    % Collect the flip time and rng values here.
    flipTime(:,end+1) = flipTimeTemp;
    rngVal(end+1) = rngValTemp;
    
    % Report what happened
    fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(correct));
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, experimentParams.nTest), correct);
    
    % Get current threshold estimate.
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %g, stderr: %g \n', 10 ^ threshold, stderr);
end

% Get flip time difference here.
flipTimeInterval = diff(flipTime);

% Close projector.
if (or(strcmp(experimentParams.runningMode,'PTB-sequential'),strcmp(experimentParams.runningMode,'PTB-directional')))
    CloseScreen;
end

%% Show results.
fprintf('%d trials recorded \n', estimator.nTrial);

% Estimate threshold and plot/report results.  This does a maximum
% likelihood based on the trials run, and is not subject to the
% discretization used by QUEST+.
%
% Return threshold value. For the mQUESTPlus Weibull PFs, the first
% parameter of the PF fit is the 0.81606 proportion correct threshold,
% when lapse rate is 0 and guess rate is 0.5. Better to make this an
% explicit parameter, however.
figure; clf;
thresholdCriterion = 0.81606;
plotSize = 10;
[threshold, para, dataOut] = estimator.thresholdMLE('showPlot', true, 'pointSize', plotSize, ...
    'thresholdCriterion', thresholdCriterion, 'returnData', true);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));
fprintf('Threshold (criterion proportion correct %0.4f): %0.2f (log10 units) / %0.4f (linear units)\n', ...
    thresholdCriterion,threshold,10^threshold);

%% Save the results.
if (SAVETHERESULTS)
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,sprintf('RunExpResults_%s_%d_cpd_%s',...
            conditionName,spatialTemporalParams.sineFreqCyclesPerDeg,dayTimestr));
        save(testFilename,'estimator','flipTime','flipTimeInterval','rngVal');
    end
end
