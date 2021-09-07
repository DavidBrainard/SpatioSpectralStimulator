%% Talk to Dropbox
% This code aims to talk to Dropbox to control other laptop from Linux box

% Initialize
clear all; close all; clc;

%% Write and save a file in Dropbox folder

filename = 'test'
savepath = '/home/colorlab/Dropbox (Aguirre-Brainard Lab)/SACC_materials/Communications/'
file = [savepath filename];
data = [1 2 3 4];

save(file,'data','-ascii');

%% Read the file

filename = 'test'
savepath = '/home/colorlab/Dropbox (Aguirre-Brainard Lab)/SACC_materials/Communications/'
file = [savepath filename];

readfile = load(file);

%% Execute a command

if readfile(1) == 1;
   
else
    
end