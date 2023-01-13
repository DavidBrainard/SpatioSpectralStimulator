%% SACC_FitCSF
%
% This is to fit CSF curve for SACC project.
%
% See also:
%    asymmetricParabolicFunc

% History:
%    1/13/23   smo    - Started on it.

%% Initialize.
clear; close all;

%% Load the data.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFilename = fullfile(testFiledir,'CSFAnalysisOutput');
    theData = load(testFilename);
    
    % Close the plots if popping up any.
    close all;
else
    error('Cannot find the data file!');
end

%% Read out the data.
%
% Subject.
subjectName = theData.subjectNameOptions;

% Threshold.
% 
% Data is aligned in (subject, SF, filter).
thresholds = theData.thresholdFitted;

% Error / Confidence interval.



%% Set variables.
myCSVals = [];
mySFVals = [3 6 9 12 18];
myWs = 1/[error or CIs];

%% Optimize the parameter p using fmincon.
myObjs = norm(myWs .* (myCSVals - asymmetricParabolicFunc(p, mySFVals)));
p_lowerBound = [10 0.5 0.1 0.1];
p_higherBound = [500 18 10 10];

p = fmincon(myObjs, p_lowerBound, p_higherBound);

%% Fitting.


%% Plot it.


%% Save the results.
