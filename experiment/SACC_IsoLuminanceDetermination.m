% SACC_IsoLuminanceDetermination
%
% This is for iso-luminance determination for SACC project. Here we display
% simple flicker (red-green) stimuli for observers to adjust the intensity
% of red until the flicker disappears.

% History:
%    08/11/22   smo     - Started on it.
%    08/17/22   smo     - Added a control part to update the intensity of
%                         red light while making red-green flicker.

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
whichChannelPrimary1 = 14;
whichChannelPrimary2 = 7;
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
frequecnyFlicker = 20; 
framesPerStim = round((1/frequecnyFlicker)/ifi); 

%% Set the image settings.
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

numButtonUp = 4;
numButtonDown = 2;
numButtonRight = 3;

primaryControlInterval = 1;

%% Start the flicker loop here.
framecounter = 0; 
while 1
    
    % End the session if the right button was pressed.
    if (stateButtonRight == false)
        stateButtonRight = Gamepad('GetButton', gamepadIndex, numButtonRight);
        if (stateButtonRight == true)
            fprintf('End the session... \n');
            break;
        end
    end
    
    % Get a gamepad response here.
    stateButtonUp = Gamepad('GetButton', gamepadIndex, numButtonUp);
    stateButtonDown = Gamepad('GetButton', gamepadIndex, numButtonDown);
    
    % Update the intensity of red light based on the key press above.
    if (stateButtonUp == true)
        % Increase the intensity of red light.
        if (intensityPrimary1 < nInputLevels-1)
            intensityPrimary1 = intensityPrimary1 + primaryControlInterval;
        end 
        fprintf('Button pressed: UP \n');
        
        % Set the button state to the initial.
        stateButtonup = false;
        
    elseif (stateButtonDown == true)
        % Decrease the intensity of red light.
        if (intensityPrimary1 > 0)
            intensityPrimary1 = intensityPrimary1 - primaryControlInterval;
        end
        fprintf('Button pressed: Down \n');
        
        % Set the button state to the initial.
        stateButtonDown = false;
    end
    
    % Update the intensity of the red light here.
    primarySetting1 = [intensityPrimary1 0 0]';
    fillColors(:,1) = primarySetting1;
    
    % Update the fill color at desired frame time.
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

%% Close the screen.
CloseScreen;
