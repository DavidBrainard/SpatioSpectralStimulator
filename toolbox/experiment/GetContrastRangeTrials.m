function [estDomainValidation preExpDataStruct] = ...
    GetContrastRangeTrials(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect, options)
% Get contrast range based on the measured threshold using Method of
% adjustment.
%
% Syntax:
%    [estDomainValidation] = GetContrastRangeTrials(...
%     sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect)
%
% Description:
%    This measures initial contrast sensitivity (threshold) using Method of
%    adjustments. In SACC project, we want to tailor test contrast range
%    per each subject. To do so, this sets the test contrast range based on
%    initial measure of contrast sensitivity.
%
%    Here we measure two contrast threshold points, one being from the
%    highest to non-visible, and the other being the lowest to visible. The
%    test contrast range is decided based on the results.
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
%    estDomainValidation                - Contrast range found based on the
%                                         trials. The number of contrast
%                                         points can be decided.
%    preExpDataStruct                   - Struct that contains all raw data
%                                         from the trials and computation
%                                         process to get the contrast range
%                                         which is saved in
%                                         estDomainValiation.
%
% Optional key/value pairs:
%    options.nContrastPoints            - Default to 6. Number of contrast
%                                         points to make for the
%                                         experiment.
%    options.higherLimThresholdEstLog   - Default to 0.3. This decides the
%                                         higher limit of the contrast
%                                         range. This number will be added
%                                         to the threshold found in log
%                                         unit.
%    options.lowerLimThresholdEstLog    - Default to -0.5. This decides the
%                                         lower limit of the contrast
%                                         range. This number will be
%                                         substracted to the threshold
%                                         found in log unit.

% History:
%    09/26/22  smo                      - Wrote it.

%% Set parameters.
arguments
    sceneParamsStruct
    experimentParams
    autoResponseParams
    window (1,1)
    windowRect (1,4)
    options.nContrastPoints (1,1) = 6
    options.higherLimThresholdEstLog (1,1) = 0.3
    options.lowerLimThresholdEstLog (1,1) = 0.5
end

%% Load null and test images.
%
% We will use only images without phase shifts here.
whichPhaseImage = 1;
nullImage = sceneParamsStruct.predefinedRGBImages{whichPhaseImage,1};
testImages = sceneParamsStruct.predefinedRGBImages(whichPhaseImage,2:end);

predefinedContrasts = sceneParamsStruct.predefinedContrasts(2:end);
predefinedContrastsLog = log10(predefinedContrasts);

%% Set the initial screen for instruction.
imageSize = size(nullImage,2);
messageInitialImage_1stLine = 'Press any button to start';
messageInitialImage_2ndLine = 'Pre-experiment session';
initialInstructionImage = insertText(nullImage,[30 imageSize/2-40; 30 imageSize/2+40],{messageInitialImage_1stLine messageInitialImage_2ndLine},...
    'fontsize',70,'Font','FreeSansBold','BoxColor',[1 1 1],'BoxOpacity',0,'TextColor','black','AnchorPoint','LeftCenter');
initialInstructionImage = fliplr(initialInstructionImage);

%% Set the contrast image for sensitivity measure.

% Make a loop to proceed this procedure, one from the highest contrast
% and the other from the lowest contrast.
%
% Start from the image with the highest contrast. Note that the first
% image is null image with out contrast pattern, so we start from
% either lowest or highest contrast.
initialContrast = [length(testImages) 1];
nInitialContrasts = length(initialContrast);

% Make variables to collect raw data.
testContrastCollect = [];
rngValuesCollect = [];
whichDirectionCollect = [];

% Start the loop for two sessions starting either from lowest or highest
% contrast.
for cc = 1:nInitialContrasts
    
    %% Display the initial screen.
    SetScreenImage(initialInstructionImage, window, windowRect,'verbose',true);
    
    % Get any button press to proceed.
    GetGamepadResp;
    disp('Practice trial is going to be started!');
    
    %% Set the initial contrast level here.
    %
    % Contrast level starts either highest one or lowest one. Starts
    % from the lower one.
    imageContrastLevel = initialContrast(cc);
    fprintf('Starting initial contrast sensitivity measure (%d/%d) \n',cc,nInitialContrasts);
    
    while 1
        % Set the contrast level.
        testContrast = predefinedContrasts(imageContrastLevel);
        fprintf('Current test contrast is = (%.4f) \n',testContrast);
        
        % Display contrast image here.
        [correct, flipTime, rngValues, whichDirectionToDisplay] = computePerformanceSACCDisplay(nullImage, testImages{imageContrastLevel}, ...
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
        
        % Collect raw data. We will not keep the data repeated.
        if ~(testContrast == testContrastCollect(end))
            testContrastCollect(end+1) = testContrast;
            rngValuesCollect(end+1) = rngValues;
            whichDirectionCollect(end+1) = whichDirectionToDisplay;
        end
        
        % Update the state according to button press.
        if strcmp(buttonPress,'down')
            fprintf('Finishing up the session... \n');
            break;
            
        elseif strcmp(buttonPress,'right')
            % Change the contrast level for next display.
            if  (cc == 1)
                % Starting from highest contrast.
                imageContrastLevel = imageContrastLevel - 1;
            elseif (cc == 2)
                % Starting from lowest contrast.
                imageContrastLevel = imageContrastLevel + 1;    
            end
            
            % Play the sound.
            MakeBeepSound('preset',correct);
            
        elseif strcmp(buttonPress,'left')
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
    thresholdFoundRawLinear(cc) = predefinedContrasts(imageContrastLevel);
    fprintf('Contrast was found at (%.3f) \n', thresholdFoundRawLinear(cc));
    fprintf('Initial contast sensitivity measure is in progress -(%d/%d) \n', cc, nInitialContrasts);
end

%% Set the contrast range based on the results for constant stimuli method.
%
% Make an average of two contrast points as initial guess of threshold.
thresholdEstLinear = mean(thresholdFoundRawLinear);
thresholdEstLog = log10(thresholdEstLinear);

% We set the contrast range based on the above threshold estimation.
higherLimTestContrastLog = thresholdEstLog + options.higherLimThresholdEstLog;
lowerLimTestContrastLog  = thresholdEstLog - options.lowerLimThresholdEstLog;

% Make contrast range equally spaced on log space.
estDomainValidationNominalLinear = logspace(lowerLimTestContrastLog, higherLimTestContrastLog, options.nContrastPoints);
estDomainValidationNominalLog = log10(estDomainValidationNominalLinear);

% Find the neareast point for each target contrast point.
for tt = 1:length(estDomainValidationNominalLog)
    [val idx] = min(abs(estDomainValidationNominalLog(tt)-predefinedContrastsLog));
    estDomainValidation(tt) = predefinedContrastsLog(idx);
end

% Delete the same values.
estDomainValidation = unique(estDomainValidation);

% Save the raw data in struct.
preExpDataStruct.rawdata.testContrast = testContrastCollect;
preExpDataStruct.rawdata.rngValues = rngValuesCollect;
preExpDataStruct.rawdata.whichDirectionToDisplay = whichDirectionCollect;

preExpDataStruct.thresholdFoundRawLinear = thresholdFoundRawLinear;
preExpDataStruct.thresholdEstLinear = thresholdEstLinear;
preExpDataStruct.estDomainValidationNominalLinear = estDomainValidationNominalLinear;

end
