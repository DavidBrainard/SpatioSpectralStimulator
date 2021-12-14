%% Initialize 
clear all; close all; clc;
tbUse('BrainardLabBase');
tbUse('SpatioSpectralStimulator');
tbUseProject('SpatioSpectralStimulator'); % Hook to the Dropbox for saving the calibration file

addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx')); % add VPixx toolbox ('Datapixx') to path (incl. all subfolders)
