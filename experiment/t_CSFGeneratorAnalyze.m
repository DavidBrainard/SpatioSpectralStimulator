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
%    02/28/22  smo -          Started on it.

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
PF = 'weibull';
conditionName = 'LminusMSmooth';

%% Load the data and PF fitting.
%
% You can loop it if you want to fit multiple data.
nData = 10;
for dd = 1:nData
    % Load the data.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('RunExpResults_%s',conditionName),'olderDate',dd-1);
        theData = load(testFilename);
    else
        error('Cannot find data file');
    end
    
    % Set fitting params and pull out the data.
    nTrials = theData.estimator.nRepeat;
    thresholdCriterion = 0.81606;
    [threshold, para, dataOut] = theData.estimator.thresholdMLE('thresholdCriterion', thresholdCriterion, 'returnData', true);
    
    % Set the contrast levels in linear unit.
    examinedContrastsLinear = 10.^dataOut.examinedContrasts;
    
    % PF fitting here.
    [paramsFitted] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, 'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE);
end
