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

%% Initialize.
clear; close all;

%% Open the projector.
initialScreenSettings = [0 0 0]';
projectorMode = true;
[window windowRect] = OpenPlainScreen(initialScreenSettings, 'projectorMode', projectorMode);

%% Set the projector LED channel settings.
%
% Make channel settings.
nPrimaries = 3;
nChannels = 16;
nInputLevels = 256;

% Set which channel to use per each primary. We will use only two
% primaries.
%
% Peak wavelength in order of channel number. Note that it is not ascending
% order.
% [422,448,476,474,506,402,532,552,558,592,610,618,632,418,658,632]
whichChannelPrimary1 = 15;
whichChannelPrimary2 = 5;
channelIntensityPrimary1 = 1;
channelIntensityPrimary2 = 1;

channelSettings = zeros(nChannels, nPrimaries);
channelSettings(whichChannelPrimary1, 1) = channelIntensityPrimary1;
channelSettings(whichChannelPrimary2, 2) = channelIntensityPrimary2;

% Set channel setting here.
SetChannelSettings(channelSettings);

%% Set parameters of the flicker.
%
% Set projector frame rate.
frameRate = 120;
ifi = 1/frameRate;

% Set the flicker frequency .
frequecnyFlicker = 40;
framesPerStim = round((1/frequecnyFlicker)/ifi);

%% Set the image settings.
%
% Set the primary colors.
intensityPrimary1 = nInputLevels-1;
intensityPrimary2 = nInputLevels/2;

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

% Set primary index to update to make a flicker.
fillColorIndexs = [1 2];
fillColorIndex = 1;

%% Set gamepad settings.
%
% Get the gamepad index for getting response.
gamepadIndex = Gamepad('GetNumGamepads');

% Set parameters for gamepad.
stateButtonUp = false;
stateButtonDown = false;
stateButtonRight = false;
actedUp = false;
actedDown = false;

numButtonUp = 4;
numButtonDown = 2;
numButtonRight = 3;

primaryControlInterval = 1;

%% Start the flicker loop here.
frameCounter = 0;

while 1
    
    % End the session if the right button was pressed.
    if (stateButtonRight == false)
        stateButtonRight = Gamepad('GetButton', gamepadIndex, numButtonRight);
        if (stateButtonRight == true)
            fprintf('Finishing up the session... \n');
            break;
        end
    end
    
    % Get a gamepad response here.
    stateButtonUp = Gamepad('GetButton', gamepadIndex, numButtonUp);
    stateButtonDown = Gamepad('GetButton', gamepadIndex, numButtonDown);
    
    % Reset acted on state when button comes back up.
    if (actedUp && ~stateButtonUp)
        actedUp = false;
    end
    if (actedDown && ~stateButtonDown)
        actedDown = false;
    end
    
    % Update the intensity of red light based on the key press above.
    if (stateButtonUp && ~actedUp)
        % Increase the intensity of red light.
        if (intensityPrimary1 < nInputLevels-1)
            intensityPrimary1 = intensityPrimary1 + primaryControlInterval;
        end
        actedUp = true;
        fprintf('Button pressed: (UP)   / Red = (%d), Green = (%d) \n', intensityPrimary1, intensityPrimary2);
        
    elseif (stateButtonDown && ~actedDown)
        % Decrease the intensity of red light.
        if (intensityPrimary1 > 0)
            intensityPrimary1 = intensityPrimary1 - primaryControlInterval;
        end
        actedDown = true;
        fprintf('Button pressed: (DOWN) / Red = (%d), Green = (%d) \n', intensityPrimary1, intensityPrimary2);
    end
    
    % Update the intensity of the red light here.
    primarySetting1 = [intensityPrimary1 0 0]';
    fillColors(:,1) = primarySetting1;
    
    % Update the fill color at desired frame time.
    if ~mod(frameCounter, framesPerStim)
        fillColorIndex = setdiff(fillColorIndexs, fillColorIndex);
        fillColor = fillColors(:,fillColorIndex);
    end
    
    % Display here.
    Screen('FillOval', window, fillColor, ovalRect);
    
    % Make a flip.
    Screen('Flip', window);
    
    % Increase frame counter.
    frameCounter = frameCounter + 1;
end

%% Close the screen.
CloseScreen;
