function PracticeTrialsOneContrastImage(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect)
% This runs practice trials with one image with the highest contrast.
%
% Syntax:
%    PracticeTrialsOneContrastImage(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect)
%
% Description:
%    It runs practice trials before starting the main experiment in SACC
%    project. It continuosly shows the same image with the highest contrast
%    so that subject can get used to how to evaluate it for the main
%    experiment.
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
%    N/A
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
messageInitialImage_2ndLine = 'Practice trials';
initialInstructionImage = insertText(nullImage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialImage_1stLine messageInitialImage_2ndLine},...
    'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
initialInstructionImage = fliplr(initialInstructionImage);

%% Display the initial screen.
SetScreenImage(initialInstructionImage, window, windowRect,'verbose',true);

% Get any button press to proceed.
GetGamepadResp;
disp('Practice trial is going to be started!');

%% Start the session.
fprintf('Practice trial starts...\n');

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
    buttonPress = GetGamepadResp;
    displayLeft = 1;
    displayRight = 2;
    
    % Update the state according to button press.
    if strcmp(buttonPress,'down')
        fprintf('Finishing up the practice trials session... \n');
        
        % Play sound as feedback when the contrast level was decided.
        numPlaySound = 3;
        for pp = 1:numPlaySound
            MakeBeepSound('preset','correct');
        end
        break;
        
    elseif strcmp(buttonPress,'left')
        if (whichDirectionToDisplay == displayLeft)
            MakeBeepSound('preset','correct');
        elseif ~(whichDirectionToDisplay == displayLeft)
            MakeBeepSound('preset','incorrect');
        end
        
    elseif strcmp(buttonPress,'right')
        if (whichDirectionToDisplay == displayRight)
            MakeBeepSound('preset','correct');
        elseif ~(whichDirectionToDisplay == displayRight)
            MakeBeepSound('preset','incorrect');
        end
    end
end

end