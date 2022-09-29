% SACC_PsychophysicalScreeningTest
%
% This is to test eligibility of subject to participate in the experiment
% based on a brief psychophysical screening test.

% History:
%    09/29/22   smo     - Started on it.

%% Initialize.
clear; close all;

%% Load the data.
%
% For now, we use random sampled date here.
data = [0.6572933 0.01733973 0.2823567 0.3686444 0.3667386 0.6939998];

%% Calculate the interquartile range.
medianData = median(data);
ICR = prctile(data,75) - prctile(data,25);

%% Print out the results.
