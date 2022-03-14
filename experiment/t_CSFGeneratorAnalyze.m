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
nData = 6;
SUBPLOT = true;
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
    [stimVec, responseVec, structVec] = combineData(theData.estimator);
    
    % Make marker size different according to the number of trials.
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
        subplot(3,2,dd); hold on;
    end 
    [paramsFitted] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
        'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE, 'figureWindow', ~SUBPLOT, 'pointSize', pointSize);
    clear pointSize;
end
