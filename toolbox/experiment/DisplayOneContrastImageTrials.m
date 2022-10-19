function DisplayOneContrastImageTrials(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect)
% This displays one contrast image with the highest contrast continuosly.
%
% Syntax:
%    DisplayOneContrastImageTrials(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect)
%
% Description:
%    This function displays one image that has the highest contrast within
%    the input images. It is written to let subject have a better idea of
%    how image is going to be displayed during the experiment.
%
% Inputs:
%    sceneParamsStruct                  - Struct containing properties of
%                                         the scene understood by this
%                                         function. It contains the RGB
%                                         null and test contrast images.
%    experimentParams                   - Struct containing properties of
%                                         running the experiments.
%    autoResponseParams                 - Parameters to get an auto
%                                         response.
%    window                             - PTB window for opened screen.
%    windowRect                         - Rect corresonding to window.
%
% Output:
%    N/A
%
% Optional key/value pairs:
%    N/A
%
% See also:
%    GetContrastRangeTrials.
%

% History:
%    10/18/22  smo                      - Wrote it.

%% Set parameters.
arguments
    sceneParamsStruct
    experimentParams
    autoResponseParams
    window (1,1)
    windowRect (1,4)
end

%% Load null and test image with highest contrast.
%
% We will use only images without phase shifts here.
whichPhaseImage = 1;
nullImage = sceneParamsStruct.predefinedRGBImages{whichPhaseImage,1};
testImage = sceneParamsStruct.predefinedRGBImages(whichPhaseImage,end);
testContrast = sceneParamsStruct.predefinedContrasts(end);

%% Set the initial screen for instruction.
imageSize = size(nullImage,2);
messageInitialImage_1stLine = 'Press any button to start';
messageInitialImage_2ndLine = 'Image display demo';
initialInstructionImage = insertText(nullImage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialImage_1stLine messageInitialImage_2ndLine},...
    'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
initialInstructionImage = fliplr(initialInstructionImage);

%% Display the initial screen.
SetScreenImage(initialInstructionImage, window, windowRect,'verbose',true);

% Get any button press to proceed.
GetGamepadResp;
disp('Practice trial is going to be started!');

%% Start the session.
fprintf('Image display demo starts...\n');

while 1
    % Set the contrast level.
    fprintf('Current test contrast is = (%.4f) \n',testContrast);
    
    % Display contrast image here.
    [correct, flipTime, rngValues, whichDirectionToDisplay] = computePerformanceSACCDisplay(nullImage, testImage{1}, ...
        sceneParamsStruct.predefinedTemporalSupport,testContrast,window,windowRect,...
        'runningMode',experimentParams.runningMode,'autoResponse',autoResponseParams,...
        'expKeyType',experimentParams.expKeyType,'beepSound',false,...
        'debugMode',experimentParams.debugMode,'movieStimuli',experimentParams.movieStimuli,...
        'movieImageDelaySec',experimentParams.movieImageDelaySec,...
        'preStimuliDelaySec',experimentParams.preStimuliDelaySec, 'addNoiseToImage', sceneParamsStruct.addNoiseToImage, ...
        'addFixationPointImage', sceneParamsStruct.addFixationPointImage,...
        'rotateImageDeg',sceneParamsStruct.rotateImageDeg, 'verbose',false);
    
    % Get a button press here.
    while 1
        buttonPress = GetGamepadResp;
        
        buttonOptions = {'right', 'down'};
        if any(strcmp(buttonPress,buttonOptions))
            break;
        end
        
        fprintf('Available button options: \n');
        fprintf('\t \t [%s] \n',buttonOptions{:});
    end
    
    % Update the state according to button press.
    if strcmp(buttonPress,'down')
        fprintf('Finishing up the demo session... \n');
       
        % Play sound as feedback when the contrast level was decided.
        numPlaySound = 3;
        for pp = 1:numPlaySound
            MakeBeepSound('preset','correct');
        end
        break;
        
    elseif strcmp(buttonPress,'right')
        % Play the sound.
        MakeBeepSound('preset','correct');
    end  
end
end
