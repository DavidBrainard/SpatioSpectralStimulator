% SACC_IsoLuminanceDetermination
%
% This is for iso-luminance determination for SACC project. Here we display
% simple flicker (red-green) stimuli for observers to adjust the intensity
% of red until the flicker disappears.

% History:
%    08/11/22   smo        - Started on it.
%    08/17/22   smo        - Added a control part to update the intensity of
%                            red light while making red-green flicker.
%    08/18/22   dhb, smo   - Now code is working fine.
%    08/19/22   smo        - Added Gaussian noise on the white background.
%    09/14/22   dhb, smo   - Flicker frequency (frame numbers) has been
%                            corrected.
%    11/08/22   smo        - Flicker session has been substituted with the
%                            function.

%% Initialize.
clear; close all;

%% Set variables.
subjectName = '002';
leftButton = 'true';

%% 1) Look at stimulus demo.
%
% For practice, subject will see the flicker stimulus by controlling the
% red starting either from top or bottom.
nTrials = 2;
GetMatchingRedForRGFlicker('nTrials',nTrials,'leftButton',leftButton);

%% 2) Practice trials.
%
% We will repeat the whole session twice so that subject will do a total of
% 4 flicker sessions (2 repeatitions x 2 starting points either top or bottom).
nTrials = 4;
dataPractice = GetMatchingRedForRGFlicker('nTrials',nTrials,'leftButton',leftButton);

%% 3) Main session.
%
% By protocol, we will meausre a total of 6 sessions (3 repeatitions x 2
% starting points either top or bottom).
nTrials = 6;
dataMain = GetMatchingRedForRGFlicker('nTrials',nTrials,'leftButton',leftButton);

%% Collect all data and save it.
theData.practice = dataPractice;
theData.main = dataMain;
