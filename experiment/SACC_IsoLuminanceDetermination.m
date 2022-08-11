% SACC_IsoLuminanceDetermination
%
% This is for iso-luminance determination for SACC project. Here we display
% simple flicker (red-green) stimuli for observers to adjust the intensity
% of red until the flicker disappears.

% History:
%    08/11/22   smo     - Started on it.

%% Initialize.
clear; close all;

%% Open the projector.
initialScreenSettings = [0 0 0]';
projectorMode = true;
[window windowRect] = OpenPlainScreen(initialScreenSettings, 'projectorMode', projectorMode);

%% Set LED channel settings here.
%
% Make channel settings.
nPrimaries = 3;
nChannels = 16;

% Set which channel to use per each primary. We will use only two
% primaries.
whichChannelPrimary1 = 14;
whichChannelPrimary2 = 7;
channelIntensityPrimary1 = 1;
channelIntensityPrimary2 = 1;

channelSettings = zeros(nChannels, nPrimaries);
channelSettings(whichChannelPrimary1, 1) = channelIntensityPrimary1;
channelSettings(whichChannelPrimary2, 2) = channelIntensityPrimary2;

% Set channel setting here.
SetChannelSettings(channelSettings);

%% Set parameters of the flicker session.
%
% Set maximum time for all stimuli to be presented in seconds.
durationSec = 5;  

% Set Hz for stimulus flicker
frequecnyFlicker = 20; 

Screen('Flip', window);
frameRate = 120;
ifi = 1/frameRate;

% Number of frames for all stimuli
framesPerFull = round(durationSec/ifi); 

% Number of frames for each stimulus
framesPerStim = round((1/frequecnyFlicker)/ifi); 

%% Make a flicker.
Framecounter = 0; %Frame counter begins at 0

% Start the flicker loop here.
while 1

    % End session
    if Framecounter==framesPerFull
        break; 
    end

    % Change background stimulus colour
    if ~mod(Framecounter,framesPerStim)
        randomcolour = rand(1, 3)*255; 
    end

    centerScreen = windowRect/2;
    centerScreenHorz = centerScreen(3);
    centerScreenVert = centerScreen(4);
    
    % ovalRect = [left top right bottom]';
    ovalSize = 100;
    ovalRect = [centerScreenHorz-ovalSize/2 centerScreenVert-ovalSize/2 ...
        centerScreenHorz+ovalSize/2 centerScreenVert+ovalSize/2]';
    Screen('FillOval', window, randomcolour, ovalRect);
    Screen('Flip', window);

    % Increase frame counter
    Framecounter = Framecounter + 1; 
end

%% Close the screen after the end of the session.
CloseScreen;
