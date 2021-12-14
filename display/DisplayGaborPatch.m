%% Display a Gabor patch on the DLP using the Psychtoolbox
% This code is to display Gabor patch for SACC project

% It is useful when your computer cannot configure the target display (LED
% screen in SACC project)

% This code can be used regardless of the configuration state as it uses
% Psychtoolbox

%% Initialize
sca;
close all;
clear all;

%% Load the Gabor patch image
% Read the gaborPatch generated from 'SpectralTestcal.m'
gaborPatch = load('testImageData1.mat'); % Load the Gabor patch image
gaborPatchImage = gaborPatch.screenSettingsImage; % Read the image data from the struct
gaborPatchImage = im2uint8(gaborPatchImage); % Convert the format of the image from 'double' to 'uin8'

%% PTB pre-setup
PsychDefaultSetup(2); 

screens = Screen('Screens');    
screenNumber = max(screens);

white = WhiteIndex(screenNumber); 
black = BlackIndex(screenNumber); 

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Read the resolution of the display

%% Display a background as plain screen
% Backgroundcolor can be set, but here it is set as the same size between
% the background and the display window, so it will not be shown on the final result
backgroundColor = [255 255 255]; % Range 0-255 (8 bits)
[window, backgroundRect] = Screen('OpenWindow', screenNumber, backgroundColor, windowRect); % Size of the background = Size of the window

%% Display the gabor patch as 'texture'
gaborPatchImage_texture = Screen('MakeTexture', window, gaborPatchImage); % Set the gaborPatchImage as 'texture' in PTB
Screen('DrawTexture', window, gaborPatchImage_texture, [], backgroundRect);
Screen('Flip', window);
