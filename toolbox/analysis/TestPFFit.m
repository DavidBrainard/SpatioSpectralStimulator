% TestCSFit.
%
% This is to test the case of fitting PF that generates negative threshold
% values when bootstrapping.
%
% See also:
%    t_CSFGeneratorAnalyze, FitPFToData.

% History:
%    4/6/23    smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set parameters to fit.
% This is the same setting that we do.
nBootstraps = 100;
paramsFree = [1 1 0 1];
minSlope = 0.1;
maxSlope = 10;
nSlopes = 20;
slopeValList = 10.^linspace(log10(minSlope),log10(maxSlope),nSlopes);

%% Target fitting data points.
% Subject 013, 18 cpd, Filter D.
stimLevels =[0.0107 0.0129 0.0188 0.0227 0.0273 0.0398 0.0481 0.0700];
pCorrect = [0.8500 0.9500 1.0000 0.9500 0.9500 1.0000 1.0000 1.0000];

%% Fitting happens here.
% Bootstrapped results are saved in 'thresholdFittedBoot'.
[paramsFitted,thresholdFitted,thresholdFittedBoot,~,~,~,~,~,~,~,~,~,~] = ...
    FitPFToData(stimLevels,pCorrect,'paramsFree', paramsFree,'nBootstraps',nBootstraps,'beta',slopeValList);

 % Set the range for the x-axis on the Figure.
 xlim([-3.3 -1]);
 
%% Check the Bootstrapped values.
I = find(thresholdFittedBoot<0);
nNegVals = length(I);
fprintf('Negative values found = (%d/%d) \n', nNegVals, nBootstraps);