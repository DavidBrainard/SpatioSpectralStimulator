% Show how to use CSF generator engine to run an experiment.
%
% Syntax:
%    t_CSFGeneratorExperiment
%
% Description:
%    Demonstrates how to use our CSF generator machinery to drive a psychophysical
%    experiment.
%
% Inputs:
%    None.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   ISETBioCSFGenerator toolbox and its tutorials.
%
% History:
%   01/18/22  dhb, so   Start writing.

%% Initialization
clear; close all;

%% Parameters
%
% Set up color direction
conditionName = 'LminusMSmooth';
colorDirectionParams = SetupColorDirection(conditionName);

%% Image spatial parameters.
%
% Image will be centered in display.
sineFreqCyclesPerDeg = 1;
gaborSdDeg = 1.5;
stimulusSizeDeg = 7;

%% Display params struct
%
% Information about the display configuration that we will use to
% initialize the scene engine but also have available to pass directly to
% the routine that runs each psychophysical trial.  This should only
% include information that we would not need to do a computational observer
% analysis, but could include things such as whether we are simulating
% trials or not, and information about the response box.  It could also
% include information required to make the ISTBio scenes.
displayParamsStruct = ...;

%% Instantiate a sceneEngine
%
% First set up the scene parameters that will be needed by
% the sceSACCDisplay.

% First step is to predefine the contrasts that we will allow the
% psychophysics to work over.  This gives us a finite list of scenes
% to compute for.
stimContrastsToTest = [];

% Then we do all of the precomputation to get the ISETBio scene sequenced for each
% possible contrast, in the color direction we are working in. We store
% this info in the sceneParams structure.  

% Initialize the scene engine. The scene params will have inside of it the
% precomputed scene sequences, temporal support, and RGB values for every
% possible contrast in our list.  So all our sceSACCDisplay funciton will
% need to do is look up the precomputed information for a particular
% contrast and return it to us.
theSceneEngine = sceneEngine(@sceSACCDisplay,sceneParams);

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
estDomain  = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
slopeRange = slopeRangeLow: slopeDelta : slopeRangeHigh;

% Note the explicit setting of the PF for the questThresholdEngine.  Using
% @qpPFWeibullLog causes it all to happen in log10 units, rather than the
% default dB units.
%
% Run 'numEstimator > 1' interleaved QUEST+ objects. In this case,
% the threshold engine calculates the running standard error (SE)
% among those objects. Essentially, this is an internal estimate of
% the precision of the current threshold estimate.  Given this, the
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
%
% Choices (comment in one):
%stopCriterion = 0.025;
stopCriterion = @(threshold, se) se / abs(threshold) < 0.01;

% Set up the estimator object.
estimator = questThresholdEngine('minTrial', 2e2, 'maxTrial', 5e3, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, ...
    'numEstimator', 4, 'stopCriterion', stopCriterion, ...
    'qpPF',@qpPFWeibullLog);

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

%% Threshold estimation with QUEST+
%
% There is some logic here to cache scenes and classifiers for each
% contrast tested, so that things will run faster if a contrast is repeated
% a second time.  This can use up space, so for real problems you may need
% to think about space-time tradeoffs, caching things to disk, etc.

% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];
while (nextFlag)
    
    % Convert log contrast -> contrast
    testContrast = 10 ^ logContrast;

    % Get the scene sequence and RGB info for the desired contrast.
    % Our scene engine provides us with the RGB values we need in the
    % status report structure.  A bit of overloading of the original
    % intent, but should work just fine.
    [theTestSceneSequences, ~, testStatusReportStruct] = ...
        theSceneEngine.compute(testContrast);

    % Run the trial and get the response. This routine
    % takes the RGB image info for the trial and returns 1 if is a
    % correct trial and 0 if incorrect trial.
    predictions = computePerformanceSACCDisplay(...
            nullStatusReportStruct, testStatusReportStruct, ...
            theSceneTemporalSupportSeconds,displayParamsStruct);
    
   
    
    % Report what happened
    fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(predictions));
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nTest), predictions);
    
    % Get current threshold estimate
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %g, stderr: %g \n', 10 ^ threshold, stderr);
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Estimate threshold and plot/report results.  This
% does a maximumu likelihood based on the trials run, and is not subject to
% the discretization used by QUEST+.
if(strcmp(questMode, 'validationMode'))
    plotSize = 50;
else
    plotSize = 10;
end

% Return threshold value. For the mQUESTPlus Weibull PFs, the first
% parameter of the PF fit is the 0.81606 proportion correct threshold,
% when lapse rate is 0 and guess rate is 0.5.  Better to make this an
% explicit parameter, however.
figure();
thresholdCriterion = 0.81606;
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', plotSize, ...
    'thresholdCriterion', thresholdCriterion);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));
fprintf('Threshold (criterion proportion correct %0.4f: %0.2f (log10 units)\n', ...
    thresholdCriterion,threshold);

%% Validation by computing the entire psychometric curve
%
% Takes a rather long time to run, change runValidation to 'true' if you
% want to run the validation
runValidation = false;
if (runValidation)
    
    logContrast = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
    pCorrect = zeros(1, length(logContrast));
    parfor idx = 1:length(logContrast)
        testContrast = 10 ^ logContrast(idx);
        [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
        
        pCorrect(idx) = mean(computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequence, ...
            theSceneTemporalSupportSeconds, 512, 512, ...
            theNeuralEngine, theRawClassifierEngine), trainFlag, testFlag);
    end
    
    hold on;
    plot(logContrast, pCorrect, '-ok', 'LineWidth', 1);
    ylim([0, 1]);
    
end

