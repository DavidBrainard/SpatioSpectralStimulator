% t_CSFGeneratorAnalyze.
%
% Syntax:
%    t_CSFGeneratorAnalyze
%
% Description:
%    This is for fitting PF to the data acquired from the experiments.
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
%    t_CSFGeneratorExperiment.

% History:
%    02/28/22  smo            Started on it.
%    03/14/22  smo            Added a plotting option to make different
%                             marker size over different number of trials.

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
PF = 'weibull';
conditionName = 'LminusMSmooth';

%% Load the data and PF fitting.
%
% You can loop it if you want to fit multiple data. nData decides the
% number of data to load from the most recent one.
%
% Set startData to 0 if you want to read the data from the most recent.
startData = 0;
nData = 1;
SUBPLOT = true;
for dd = 1:nData
    % Load the data.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('RunExpResults_%s',conditionName),'olderDate',startData+dd-1);
        theData = load(testFilename);
    else
        error('Cannot find data file');
    end
    
    % Set fitting params and pull out the data.
    nTrials = theData.estimator.nRepeat;
    [stimVec, responseVec, structVec] = combineData(theData.estimator);
    
    % Set marker size here. We will use this number to plot the results to
    % have different marker size according to the number of trials. Here
    % we used the same method to decide the size of each marker as
    % 'thresholdMLE' does.
    stimVal = unique(stimVec);
    pCorrect = zeros(1,length(stimVal));
    for idx = 1:length(stimVal)
        prop = responseVec(stimVec == stimVal(idx));
        pCorrect(idx) = sum(prop) / length(prop);
        pointSize(idx) = 10 * 100 / length(stimVec) * length(prop);
    end
    
    thresholdCriterion = 0.81606;
    [threshold, para, dataOut] = theData.estimator.thresholdMLE(...
        'thresholdCriterion', thresholdCriterion, 'returnData', true);
    
    % Set the contrast levels in linear unit.
    examinedContrastsLinear = 10.^dataOut.examinedContrasts;
    
    % PF fitting here.
    if (SUBPLOT)
        subplot(1,1,dd); hold on;
    end
    [paramsFitted(:,dd)] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
        'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE, 'figureWindow', ~SUBPLOT, 'pointSize', pointSize);
    clear pointSize;
end

%% Add the mean threshold from Adaptive method if you want to compare it with the result from Validation method.
%
% This is temporary, we may make this part more elaborate later on.
ADDMEANTHRESHOLD = false;

if (ADDMEANTHRESHOLD)
    meanThreshold = 0.0030;
    stdThreshold = 0.0010;
    plot(meanThreshold, thresholdCriterion, 'ko','MarkerFaceColor','b','MarkerSize',12);
    e = errorbar(meanThreshold,thresholdCriterion,stdThreshold,'horizontal');
    e.Color = 'blue';
end