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
nInputLevels = 256;

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
durationSec = 3;  

% Set Hz for stimulus flicker
frequecnyFlicker = 25; 
frameRate = 120;
ifi = 1/frameRate;

% Number of frames for all stimuli
framesPerFull = round(durationSec/ifi); 

% Number of frames for each stimulus
framesPerStim = round((1/frequecnyFlicker)/ifi); 

%% Make a flicker.
%
% Set the primary colors.
intensityPrimary1 = nInputLevels-1;
intensityPrimary2 = nInputLevels-1;

primarySetting1 = [intensityPrimary1 0 0]';
primarySetting2 = [0 intensityPrimary2 0]';
fillColors = [primarySetting1 primarySetting2];

% Set the position of the circle on the screen.
% Note that ovalRect = [left top right bottom]';
centerScreen = windowRect/2;
centerScreenHorz = centerScreen(3);
centerScreenVert = centerScreen(4);

ovalSize = 150;
ovalRect = [centerScreenHorz-ovalSize/2 centerScreenVert-ovalSize/2 ...
    centerScreenHorz+ovalSize/2 centerScreenVert+ovalSize/2]';

% Frame counter begins at 0
framecounter = 0; 

% Set primary index to update to make a flicker.
fillColorIndexs = [1 2];
fillColorIndex = 1;

% Start the flicker loop here.
while 1
    
    % End session
    if framecounter == framesPerFull;
%         numButtonRight = 3;
%         responseGamePad = GetGamepadResp2AFC('numButtonB',numButtonRight,'verbose',true); 
        break;
    end
    
    % Get a response. Here we change the intensity of red light.
%     if
%         numButtonUp    = 4;
%         numButtonDown  = 2;
%         pressButtonA = 1;
%         pressButtonB = 2;
%         responseGamePad = GetGamepadResp2AFC('numButtonA', numButtonUp, 'numButtonB',numButtonDown,'verbose',true);
%         
%         switch responseGamePad
%             case pressButtonA
%                 % Increase the intensity of red light.
%                 intensityPrimary1 = intensityPrimary1 + 1;
%             case pressButtonB
%                 % Decrease the intensity of red light.
%                 intensityPrimary1 = intensityPrimary1 - 1;
%         end
%         primarySetting1 = [intensityPrimary1 0 0]';
%         fillColors = [primarySetting1 primarySetting2];
%     end
    
     if ~mod(framecounter, framesPerStim)
        fillColorIndex = setdiff(fillColorIndexs, fillColorIndex);
        fillColor = fillColors(:,fillColorIndex); 
     end
    
    % Display here.
    Screen('FillOval', window, fillColor, ovalRect);
    Screen('Flip', window);

    % Increase frame counter.
    framecounter = framecounter + 1; 
end

%% Close the screen after the end of the session.
CloseScreen;
