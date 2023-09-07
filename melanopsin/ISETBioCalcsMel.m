% ISETBioCalcsMel

% Initialize
clear; close all;

% Load the data
testFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
testFilename = fullfile(testFiledir,'testImageDataISETBio');
theData = load(testFilename);
disp('Data loaded');

