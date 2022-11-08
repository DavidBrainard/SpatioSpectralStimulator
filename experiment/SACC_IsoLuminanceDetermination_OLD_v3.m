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
    
    framesPerStimIndexs = [1 2 1 2];
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
plainImageRed = plainImageBase;
intensityPrimary1 = 180;
plainImageRed(:,:,1) = intensityPrimary1;

% Green plain image.
plainImageGreen = plainImageBase;
intensityPrimary2 = 75;
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
actedRight = false;
actedLeft = false;

numButtonUp = 4;
numButtonDown = 2;
numButtonRight = 3;
numButtonLeft = 1;

primaryControlInterval = 4;

%% Make PTB texture for all possible settings.
%
% As making texture takes an extra time which could cause delay for
% diplaying a texture on desired time, so here we make all textures before
% starting the loop.
%
% Red.
fillColorTemp = plainImageBase;
BACKGROUND = 'white';

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
    
    [imageTextureRed(pp), imageWindowRect] = MakeImageTexture(fillColorTemp, window, windowRect,...
        'addFixationPointImage','crossbar','verbose',false);
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
imageTextureGreen = MakeImageTexture(fillColorGreen, window, windowRect,...
    'addFixationPointImage','crossbar','verbose',false);

%% Make a loop here.
redStartingPointOptions = {'top' 'bottom'};
idxRedStartingPoint = [1 2 1 2];
nTrials = length(idxRedStartingPoint);

for tt = 1:nTrials
    redStartingPoint = redStartingPointOptions{idxRedStartingPoint(tt)};
    
    switch redStartingPoint
        case 'top'
            intensityPrimary1 = 180;
        case 'bottom'
            intensityPrimary1 = 55;
    end
    
    %% Start the flicker loop here.
    frameCounter = 1;
    nextFrame = 1;
    frameIndexCounter = 1;
    
    % Back to the default setting.
    stateButtonRight = false;
    stateButtonDown = false;
    
    %% Set the initial screen for instruction.
    plainWhiteImage = plainImageBase;
    intensityWhite = 1;
    plainWhiteImage(:,:,:) = intensityWhite;
    
    imageSize = stimulusN;
    messageInitialImage_1stLine = 'Press any button to start';
    messageInitialImage_2ndLine = append('Flicker session (',num2str(tt),'/',num2str(nTrials),')');
    initialInstructionImage = insertText(plainWhiteImage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialImage_1stLine messageInitialImage_2ndLine},...
        'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
    initialInstructionImage = fliplr(initialInstructionImage);
    
    % Display the initial screen.
    SetScreenImage(initialInstructionImage, window, windowRect,'verbose',true);
    
    % Get any button press to proceed.
    GetGamepadResp;
    disp('Flicker session is going to be started!');
    
    %% Start a flicker loop here.
    while 1
        
        % End the session if the right button was pressed.
        if (stateButtonDown == false)
            stateButtonDown = Gamepad('GetButton', gamepadIndex, numButtonDown);
            if (stateButtonDown == true)
                fprintf('Finishing up the session... \n');
                break;
            end
        end
    
        % Get a gamepad response here.
        stateButtonRight = Gamepad('GetButton', gamepadIndex, numButtonRight);
        
        % Reset acted on state when button comes back up.
        if (actedRight && ~stateButtonRight)
            actedRight = false;
        end
        
        % Update the intensity of red light based on the key press above.
        if (stateButtonRight && ~actedRight)
            if strcmp(redStartingPoint,'top')
                % Decrease the intensity of red light.
                if (intensityPrimary1 > 0)
                    intensityPrimary1 = intensityPrimary1 - primaryControlInterval;
                end
                % Cut the value on the range.
                if (intensityPrimary1 < 0)
                    intensityPrimary1 = 0;
                end
                
            elseif strcmp(redStartingPoint,'bottom')
                % Increase the intensity of red light.
                if (intensityPrimary1 < nInputLevels-1)
                    intensityPrimary1 = intensityPrimary1 + primaryControlInterval;
                end
                % Cut the value over the maximum.
                if (intensityPrimary1 > nInputLevels-1)
                    intensityPrimary1 = nInputLevels-1;
                end
            end
            actedRight = true;
            fprintf('Button pressed! Red = (%d), Green = (%d) \n', intensityPrimary1, intensityPrimary2);
        end
        
        % Update the intensity of the red light here.
        imageTextures = [imageTextureRed(intensityPrimary1+1) imageTextureGreen];
        
        % Update the fill color at desired frame time.
        if (frameCounter >= nextFrame)
            fillColorIndex = setdiff(fillColorIndexs, fillColorIndex);
            imageTexture = imageTextures(fillColorIndex);
            
            % Change the frames per stim if there are more than one target frames.
            if ~isempty(framesPerStimSet)
                framesPerStimIndex = framesPerStimIndexs(frameIndexCounter);
                framesPerStim = framesPerStimSet(framesPerStimIndex);
                
                % Update the frame index counter here.
                if (frameIndexCounter < length(framesPerStimIndexs))
                    frameIndexCounter = frameIndexCounter + 1;
                else
                    % Set the counter back to 1 if it ran one set of cycle.
                    frameIndexCounter = 1;
                end
            end
            nextFrame = frameCounter+framesPerStim;
        end
        
        % Make a flip.
        flipTime(frameCounter) = FlipImageTexture(imageTexture, window, imageWindowRect, 'verbose', false);
        
        % Count the frame.
        frameCounter = frameCounter + 1;
    end
    
    % Finishing up with the beep sounds.
    nBeeps = 3;
    for bb = 1:nBeeps
        MakeBeepSound('preset','correct');
    end
    
    % Print out the matching results.
    fprintf('The matching intensity of red = (%d) \n', intensityPrimary1);
    
    % Save the matching value.
    data(tt) = intensityPrimary1;
    fprintf('The matching result has been saved! \n');
end

%% Save the results as a file.


%% Close the screen.
CloseScreen;
