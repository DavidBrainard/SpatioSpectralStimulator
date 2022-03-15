% t_CSFMergeDataAndFitPF
%
% This is to merge the data more than two conditions and fit PF to the
% data.
%
% Syntax:
%    t_CSFMergeDataAndFitPF
%
% Description:
%    This is to merge the data more than two conditions and fit PF to the
%    data.
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
%    t_CSFGeneratorAnalyze.

% History:
%    03/15/22  smo            Started on it.

%% Initialize.
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
% 18 cpd, adaptive; startData = 1  / nData = 5
% 3  cpd, adaptive; startData = 7  / nData = 5
% 1  cpd, adpative; startData = 12 / nData = 5
startData = 1;
nData = 5;
SUBPLOT = false;
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
    [stimVec(:,dd), responseVec(:,dd), structVec(:,dd)] = combineData(theData.estimator);
end

% Reshape the data.
stimVec_reshape = reshape(stimVec,[1 size(stimVec,1)*size(stimVec,2)]);
responseVec_reshape = reshape(responseVec,[1 size(responseVec,1)*size(responseVec,2)]);

% Set marker size here. We will use this number to plot the results to
% have different marker size according to the number of trials. Here
% we used the same method to decide the size of each marker as
% 'thresholdMLE' does.
stimVal = unique(stimVec_reshape);
pCorrect = zeros(1,length(stimVal));
for idx = 1:length(stimVal)
    prop = responseVec_reshape(stimVec_reshape == stimVal(idx));
    pCorrect(idx) = sum(prop) / length(prop);
    pointSize(idx) = 10 * 100 / length(stimVec) * length(prop);
end

thresholdCriterion = 0.81606;
[threshold, para, dataOut] = theData.estimator.thresholdMLE(...
    'thresholdCriterion', thresholdCriterion, 'returnData', true);

% Set the contrast levels in linear unit.
examinedContrastsLinear = 10.^stimVal;

% PF fitting here.
[paramsFitted] = FitPFToData(examinedContrastsLinear, pCorrect, ...
    'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE, 'figureWindow', ~SUBPLOT, 'pointSize', pointSize);
