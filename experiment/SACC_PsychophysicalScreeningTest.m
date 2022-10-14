% SACC_PsychophysicalScreeningTest
%
% This is to test eligibility of subject to participate in the experiment
% based on a brief psychophysical screening test.

% History:
%    09/29/22   smo     - Started on it.

%% Initialize.
clear; close all;

%% Get the current dir.
dir = cd;

%% Load the data.
%
% Set the dir.
userName = 'Mela Nopsin';
expName = 'Exp_0001J';
subjectID = 'SACC_002';
fileDir = sprintf('~/Aguirre-Brainard Lab Dropbox/%s/MTRP_data/%s/Subject_%s',userName,expName,subjectID);
cd(fileDir);

% Load the file.
fileName1 = sprintf('%s_1',subjectID);
fileName2 = sprintf('%s_2',subjectID);
theData1 = readtable(fileName1);
theData2 = readtable(fileName2);

% Export the threshold data.
numVarThreshold = 3;
thresholds1 = theData1{:,numVarThreshold};
thresholds2 = theData2{:,numVarThreshold};
thresholds = [thresholds1; thresholds2];

%% Test proficiency here.
%
% Subjects will be excluded for a median threshold of > 1.5% or an
% interquartile range across the 12 measures of > 2.0%. 
% Record as pass/fail.
criteriaMedianThreshold = 1.5;
criteriaICR = 2;

medianThreshold = median(thresholds);
ICR = prctile(thresholds,75) - prctile(thresholds,25);

% Print out the results.
if or(medianThreshold > criteriaMedianThreshold, ICR > criteriaICR)
    fprintf('    FAIL! Median threshold = %.3f / ICR = %.3f \n', medianThreshold, ICR);
else
    fprintf('    PASS! Median threshold = %.3f / ICR = %.3f \n', medianThreshold, ICR);
end

%% Back to current dir.
cd(dir);
