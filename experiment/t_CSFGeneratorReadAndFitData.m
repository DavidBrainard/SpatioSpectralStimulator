% t_CSFGeneratorReadAndFitData.
%
% Syntax:
%    t_CSFGeneratorReadAndFitData
%
% Description:
%    This is to read and PF fit the data for CSF generator experiment.
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
%    t_CSFGeneratorExperiment, t_CSFGeneratorAnalyze

% History:
%    08/30/22  smo    - Wrote it.

%% Start over.
clear; close all;

%% Set parameters here.
%
% Note that the most recent experiment with David was on 2022-08-26.
%
% Change 'sineFreqCyclesPerDeg' within [3, 6, 9, 12, 18] to see the
% results over different spatial frequency.
VERBOSE = true;
PF = 'weibull';
conditionName = 'LminusMSmooth';
sineFreqCyclesPerDeg = 9;
fileDate = '2022-08-26';

%% Load the data.
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = GetMostRecentFileName(testFiledir,...
        sprintf('RunExpResults_%s_%d_cpd_%s',conditionName,sineFreqCyclesPerDeg,fileDate),'olderDate',0);
    theData = load(testFilename);
else
    error('Cannot find data file');
end

% Pull out the data here.
nTrials = theData.estimator.nRepeat;
[stimVec, responseVec, structVec] = combineData(theData.estimator);

%% PF fitting and print out the results.
thresholdCriterion = 0.81606;
[threshold, para, dataOut] = theData.estimator.thresholdMLE(...
    'thresholdCriterion', thresholdCriterion, 'returnData', true);

% Set the contrast levels in linear unit.
examinedContrastsLinear = 10.^dataOut.examinedContrasts;

% PF fitting happens here.
[paramsFitted] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
    'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE);
