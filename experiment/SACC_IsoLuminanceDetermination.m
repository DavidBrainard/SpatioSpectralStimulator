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
whichChannelPrimary3 = [1:16];

channelIntensityPrimary1 = 1;
channelIntensityPrimary2 = 1;
channelIntensityPrimary3 = 0;

channelSettings = zeros(nChannels, nPrimaries);
channelSettings(whichChannelPrimary1, 1) = channelIntensityPrimary1;
channelSettings(whichChannelPrimary2, 2) = channelIntensityPrimary2;
channelSettings(whichChannelPrimary3, 3) = channelIntensityPrimary3;

% Set channel setting here.
SetChannelSettings(channelSettings);

%% Set parameters of the flicker.
%
% Set projector frame rate.
frameRate = 120;
ifi = 1/frameRate;

% Set the flicker frequency .
frequecnyFlicker = 25;
framesPerStim = round((1/frequecnyFlicker)/ifi);
framesPerStim = framesPerStim/2;

% If the frame number is not integer, take two closest from the number.
if ~(framesPerStim == ceil(framesPerStim))
    framesPerStimSet = [floor(framesPerStim) ceil(framesPerStim)];
    framesPerStim = framesPerStimSet(1);
    
    framesPerStimIndexs = [1 2 2 1];
    framesPerStimIndex = framesPerStimIndexs(1);
else
    framesPerStimSet = [];
end

%% Make Gaussian window and normalize its max to one.
%
% Set the window parameters. ScreenPixelsPerDeg is based on the Maxwellian
% view settings for SACC project.
stimulusN = 996;
centerN = stimulusN/2;
gaborSdDeg = 0.5;
screenPixelsPerDeg = 142.1230;
gaborSdPixels = gaborSdDeg * screenPixelsPerDeg;

% Make a gaussian window here.
gaussianWindowBase = zeros(stimulusN, stimulusN, 3);
gaussianWindow = normpdf(MakeRadiusMat(stimulusN,stimulusN,centerN,centerN),0,gaborSdPixels);
gaussianWindowBGBlack = gaussianWindow/max(gaussianWindow(:));
gaussianWindowBGWhite = 1 - gaussianWindowBGBlack;

% Make plain red/green images here.
plainImageBase = zeros(stimulusN, stimulusN, 3);

% Red plain image.
redStartingPoint = 'bottom';

plainImageRed = plainImageBase;
switch redStartingPoint
    case 'top'
        intensityPrimary1 = nInputLevels-1;
    case 'bottom'
        intensityPrimary1 = 0;
end
plainImageRed(:,:,1) = intensityPrimary1;

% Green plain image.
plainImageGreen = plainImageBase;
intensityPrimary2 = round(nInputLevels * 0.3);
plainImageGreen(:,:,2) = intensityPrimary2;

% Set primary index to update to make a flicker.
fillColors = {plainImageRed plainImageGreen};
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
stateButtonLeft = false;

actedUp = false;
actedDown = false;
actedLeft = false;

numButtonUp = 4;
numButtonDown = 2;
numButtonRight = 3;
numButtonLeft = 1;

primaryControlIntervals = [1 10];
primaryControlIntervalIndexs = [1 2];
primaryControlIntervalIndex = 2;

primaryControlInterval = primaryControlIntervals(primaryControlIntervalIndex);

%% Make PTB texture for all possible settings.
%
% As making texture takes an extra time which could cause delay for
% diplaying a texture on desired time, so here we make all textures before
% starting the loop.
%
% Red.
fillColorTemp = plainImageBase;
BACKGROUND = 'black';

% Make a loop to make all possible image textures here.
for pp = 1:nInputLevels
    fillColorTemp(:,:,1) = pp-1;
    fillColorTemp = fillColorTemp./(nInputLevels-1);
    
    % Add Gaussian window here.
    switch BACKGROUND
        case 'white'
            fillColorTemp = fillColorTemp + gaussianWindowBGWhite;
        case 'black'
            fillColorTemp = fillColorTemp .* gaussianWindowBGBlack;
    end
    
    [imageTextureRed(pp), imageWindowRect] = MakeImageTexture(fillColorTemp, window, windowRect,'verbose',false);
    fprintf('Image texture has been created - (%d/%d) \n', pp, nInputLevels);
end

% Green.
plainImageGreen = plainImageGreen./(nInputLevels-1);
switch BACKGROUND
    case 'white'
        fillColorGreen = plainImageGreen + gaussianWindowBGWhite;
    case 'black'
        fillColorGreen = plainImageGreen .* gaussianWindowBGBlack;
end
imageTextureGreen = MakeImageTexture(fillColorGreen, window, windowRect,'verbose',false);

%% Start the flicker loop here.
frameCounter = 0;
frameIndexCounter = 1;

% Start a flicker loop here.
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
    stateButtonLeft = Gamepad('GetButton', gamepadIndex, numButtonLeft);
    
    % Reset acted on state when button comes back up.
    if (actedUp && ~stateButtonUp)
        actedUp = false;
    end
    if (actedDown && ~stateButtonDown)
        actedDown = false;
    end
    if (actedLeft && ~stateButtonLeft)
        actedLeft = false;
    end
    
    % Update the intensity of red light based on the key press above.
    if (stateButtonUp && ~actedUp)
        % Increase the intensity of red light.
        if (intensityPrimary1 < nInputLevels-1)
            intensityPrimary1 = intensityPrimary1 + primaryControlInterval;
        end
        % Cut the value over the maximum.
        if (intensityPrimary1 > nInputLevels-1)
            intensityPrimary1 = nInputLevels-1;
        end
        actedUp = true;
        fprintf('Button pressed: (UP)   / Red = (%d), Green = (%d) \n', intensityPrimary1, intensityPrimary2);
        
    elseif (stateButtonDown && ~actedDown)
        % Decrease the intensity of red light.
        if (intensityPrimary1 > 0)
            intensityPrimary1 = intensityPrimary1 - primaryControlInterval;
        end
       
        % Cut the value on the range.
        if (intensityPrimary1 < 0)
            intensityPrimary1 = 0;
        end
        
        actedDown = true;
        fprintf('Button pressed: (DOWN) / Red = (%d), Green = (%d) \n', intensityPrimary1, intensityPrimary2);
        
    elseif (stateButtonLeft && ~actedLeft)
        % Update the interval as desired.
        primaryControlIntervalIndex = setdiff(primaryControlIntervalIndexs, primaryControlIntervalIndex);
        primaryControlInterval = primaryControlIntervals(primaryControlIntervalIndex);
        actedLeft = true;
        fprintf('Button pressed: (LEFT) / Control interval is now = (%d) \n', primaryControlInterval);
        
        % Play a sound when the interval changes.
        switch primaryControlInterval
            case primaryControlIntervals(1)
                beepSound = 'incorrect';
            case primaryControlIntervals(2)
                beepSound = 'correct';
        end
        MakeBeepSound('preset',beepSound);
    end
   
     % Update the intensity of the red light here.
    imageTextures = [imageTextureRed(intensityPrimary1+1) imageTextureGreen];   
            
    % Update the fill color at desired frame time.
    if ~mod(frameCounter, framesPerStim)
        fillColorIndex = setdiff(fillColorIndexs, fillColorIndex);
        imageTexture = imageTextures(fillColorIndex);
        
        % Change the frames per stim if there are more than one target frames.
        if ~isempty(framesPerStimSet)
            framesPerStimIndex = framesPerStimIndexs(frameIndexCounter);
            framesPerStim = framesPerStimSet(framesPerStimIndex);
            
            % Update the frame counter index here.
            if (frameIndexCounter < length(framesPerStimIndexs))
                frameIndexCounter = frameIndexCounter + 1;
            else
                % Set the counter back to 1 if it ran one set of cycle.
                frameIndexCounter = 1;
            end
        end
    end
    
    % Make a flip.
    flipTime(frameCounter+1) = FlipImageTexture(imageTexture, window, imageWindowRect, 'verbose', false);
    
    % Count the frame.
    frameCounter = frameCounter + 1;
end

%% Close the screen.
CloseScreen;

% Print out the matching results.
fprintf('The matching intensity of red = (%d) \n', intensityPrimary1);
