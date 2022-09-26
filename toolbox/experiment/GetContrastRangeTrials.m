function [estDomainValidation] = GetContrastRangeTrials(...
    sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect, options)
% Get contrast range for Constant stimuli method using Method of
% adjustment.
%
% Syntax:
%    [estDomainValidation] = GetContrastRangeTrials(nullImage, testImages, ...
%     experimentParams, sceneParamsStruct, autoResponseParams)
%
% Description:
%    Measure initial contrast sensitivity using Method of adjustments. We
%    will use method of adjustments to make an initial measure of each
%    subject's contrast sensitivity. We will set the contrast range based
%    on this response.
%
% Inputs:
%    experimentParams
%    sceneParamsStruct
%    autoResponseParams
%    window
%    windowRect
%
% Output:
%    estDomainValidation     -
%
% Optional key/value pairs:
%    options.nContrastPoints

% History:
%    09/26/22  smo           - Wrote it.

%% Set parameters.
arguments
    sceneParamsStruct 
    experimentParams 
    autoResponseParams
    window (1,1)
    windowRect (1,4)
    options.nContrastPoints (1,1) = 6
end

%%
nullImage = sceneParamsStruct.predefinedRGBImages{1,1};
testImages = sceneParamsStruct.predefinedRGBImages(1,2:end);

%% Set the initial screen for instruction.
imageSize = size(nullImage,2);
messageInitialImage_1stLine = 'Press any button to start';
messageInitialImage_2ndLine = 'Pre-experiment session';
initialImage = insertText(nullImage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialImage_1stLine messageInitialImage_2ndLine},...
    'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
initialImage = fliplr(initialImage);

%% Set the contrast image for sensitivity measure.
initialTestImage = testImages(1,:);

% Get the gamepad index for getting response.
gamepadIndex = Gamepad('GetNumGamepads');

% Allocate the button numbers.
numButtonUp    = 4;
numButtonRight = 3;
numButtonLeft  = 1;
numButtonDown  = 2;

% Make a loop to proceed this procedure, one from the highest contrast
% and the other from the lowest contrast.
%
% Start from the image with the highest contrast. Note that the first
% image is null image with out contrast pattern, so we start from
% either 2nd or 20th (which is the highest contrast).
initialImageContrastLevels = [length(initialTestImage) 1];
nInitialImageContrastlevels = length(initialImageContrastLevels);

for cc = 1:nInitialImageContrastlevels
    
    %% Display the initial screen.
    SetScreenImage(initialImage, window, windowRect,'verbose',true);
    
    % Get a button press to proceed.
    if (strcmp(experimentParams.expKeyType,'gamepad'))
        switch (experimentParams.runningMode)
            case 'PTB-sequential'
                responseGamePad = GetGamepadResp2AFC('verbose',true);
            case 'PTB-directional'
                if (sceneParamsStruct.rotateImageDeg == 0)
                    responseGamePad  = GetGamepadResp2AFC('numButtonA', numButtonUp, 'numButtonB',numButtonRight,'verbose',true);
                else
                    responseGamePad  = GetGamepadResp2AFC('numButtonA', numButtonLeft, 'numButtonB',numButtonRight,'verbose',true);
                end
        end
        possibleResponseGamePad = [numButtonUp numButtonRight numButtonLeft];
        if (any(responseGamePad == possibleResponseGamePad))
            disp('Practice trial is going to be started!');
        end
    end
    
    %% Set the starting contrast level here.
    %
    % Contrast level starts either highest one or lowest one. Starts
    % from the lower one.
    imageContrastLevel = initialImageContrastLevels(cc);
    
    % Show starting message.
    fprintf('Starting initial contrast sensitivity measure (%d/%d) \n',cc,nInitialImageContrastlevels);
    
    while 1
        % Set the initial button press state.
        stateButtonUp = false;
        stateButtonDown = false;
        stateButtonRight = false;
        stateButtonLeft = false;
        
        % Set the contrast level.
        initialMeasureTestContrast = sceneParamsStruct.predefinedContrasts(imageContrastLevel);
        fprintf('Current test contrast is = (%.4f) \n',initialMeasureTestContrast);
        
        % Set auto response params.
        autoResponseParams.psiFunc = @qpPFWeibullLog;
        autoResponseParams.thresh = 0.004;
        autoResponseParams.slope = 2;
        autoResponseParams.guess = 0.5;
        autoResponseParams.lapse = 0.01;
        autoResponseParams.psiParams = [log10(autoResponseParams.thresh) autoResponseParams.slope autoResponseParams.guess autoResponseParams.lapse];
        
        % Display contrast image here.
        [correct] = computePerformanceSACCDisplay(nullImage, initialTestImage{imageContrastLevel}, ...
            sceneParamsStruct.predefinedTemporalSupport,sceneParamsStruct.predefinedTemporalSupportCrossbar,initialMeasureTestContrast,window,windowRect,...
            'runningMode',experimentParams.runningMode,'autoResponse',autoResponseParams,...
            'expKeyType',experimentParams.expKeyType,'beepSound',false,...
            'debugMode',experimentParams.debugMode,'movieStimuli',experimentParams.movieStimuli,...
            'movieImageDelaySec',experimentParams.movieImageDelaySec,...
            'preStimuliDelaySec',experimentParams.preStimuliDelaySec, 'addNoiseToImage', sceneParamsStruct.addNoiseToImage, ...
            'addFixationPointImage', sceneParamsStruct.addFixationPointImage,...
            'rotateImageDeg',sceneParamsStruct.rotateImageDeg, 'verbose',false);
        
        % Waiting for a button press to continue or finish the session.
        % End the session if the right button was pressed.
        while (stateButtonLeft == false && stateButtonDown == false && stateButtonRight == false)
            stateButtonLeft = Gamepad('GetButton', gamepadIndex, numButtonLeft);
            stateButtonDown = Gamepad('GetButton', gamepadIndex, numButtonDown);
            stateButtonRight = Gamepad('GetButton', gamepadIndex, numButtonRight);
        end
        
        if (stateButtonDown)
            fprintf('Finishing up the session... \n');
            break;
            
        elseif (stateButtonRight)
            % Change the contrast level for next display.
            if (cc == 2)
                imageContrastLevel = imageContrastLevel + 1;
            elseif (cc == 1)
                imageContrastLevel = imageContrastLevel - 1;
            end
            
            % Play the sound.
            MakeBeepSound('preset',correct);
            
        elseif (stateButtonLeft)
            % Show the same contrast level again for next display.
            imageContrastLevel = imageContrastLevel;
            
            % Play the feedback sound.
            numPlaySound = 2;
            for pp = 1:numPlaySound
                MakeBeepSound('preset',correct);
            end
        end
    end
    
    % Play sound as feedback when the contrast level was decided.
    numPlaySound = 3;
    for pp = 1:numPlaySound
        MakeBeepSound('preset',correct);
    end
    
    % Print out the contrast level we found.
    contrastFound(cc) = sceneParamsStruct.predefinedContrasts(imageContrastLevel);
    fprintf('Contrast was found at (%.3f) \n', contrastFound(cc));
    fprintf('Initial contast sensitivity measure has been finished!-(%d/%d) \n', cc, nInitialImageContrastlevels);
end

%% Set the contrast range based on the results for constant stimuli method.
%
% Make an average of two contrast points as initial guess of threshold.
thresholdEstLinear = mean(contrastFound);
thresholdEstLog = log10(thresholdEstLinear);

highLimitContrastLog = thresholdEstLog + 0.3;
lowLimitContrastLog  = thresholdEstLog - 0.5;

options.nContrastPoints = 6;
estDomainValidationLinear = logspace(lowLimitContrastLog, highLimitContrastLog, options.nContrastPoints);

% Set the contrast range for Constant stimuli (Validation) method.
estDomainValidationLogNominal = log10(estDomainValidationLinear);
predefinedContrastsLog = log10(sceneParamsStruct.predefinedContrasts);

% Find the nearest contrast within the predefined contrast range.
for tt = 1:length(estDomainValidationLogNominal)
    [val idx] = min(abs(estDomainValidationLogNominal(tt)-predefinedContrastsLog));
    estDomainValidation(tt) = predefinedContrastsLog(idx);
end
end
