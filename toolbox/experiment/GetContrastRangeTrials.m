function [estDomainValidation preExpDataStruct] = ...
    GetContrastRangeTrials(sceneParamsStruct, experimentParams, autoResponseParams, window, windowRect, options)
% Get contrast range based on the measured threshold using Method of
% adjustment.
%
% Syntax:
%    [estDomainValidation preExpDataStruct] = GetContrastRangeTrials(...
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
    options.replay (1,1) = false
    options.nContrastPoints (1,1) = 8
    options.higherLimThresholdEstLog (1,1) = 0.4
    options.lowerLimThresholdEstLog (1,1) = 0.4
    options.verbose (1,1) = true
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
% highest contrast.
%
% Eventually we want to get a total of four initial threshold esitmate, two
% starting from the high and the other two from the low. So, here we simply
% set the initial contrast as four number in an array that contains two
% contrasts (highest and lowest) and doubled it.
idxContrastHigh = length(testImages);
idxContrastLow = 1; 
idxInitialContrast = [idxContrastHigh idxContrastLow idxContrastHigh idxContrastLow];
nInitialContrasts = length(idxInitialContrast);

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
    imageContrastLevel = idxInitialContrast(cc);
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
        while 1
            buttonPress = GetGamepadResp;
            if (options.replay)
                if any(strcmp(buttonPress,{'left','right','down'}))
                    break
                end
            else
                if any(strcmp(buttonPress,{'right','down'}))
                    break
                end
            end
        end
        
        % Collect raw data. We will not keep the data repeated.
        if ~(strcmp(buttonPress,'left'))
            testContrastCollect(end+1) = testContrast;
            rngValuesCollect{end+1} = rngValues;
            whichDirectionCollect(end+1) = whichDirectionToDisplay;
        end
        
        % Update the state according to button press.
        if strcmp(buttonPress,'down')
            fprintf('Finishing up the session... \n');
            break;
            
        elseif strcmp(buttonPress,'right')
            
            % Change the contrast level for next display.
            %
            % Case 1) Starting from the highest contrast.
            if  (idxInitialContrast(cc) == idxContrastHigh)                
                if ~(imageContrastLevel == idxContrastLow)
                    imageContrastLevel = imageContrastLevel - 1;
                    
                    % Play the sound.
                    MakeBeepSound('preset','correct');                   
                else
                    % Print out it reached the limit.                
                    fprintf('\t You reached the lowest contrast! \n');
                    
                    % Play the sound for feedback.
                    MakeBeepSound('preset','incorrect');
                end
                
                % Case 2) Starting from the lowest contrast.
            elseif (idxInitialContrast(cc) == idxContrastLow)
                if ~(imageContrastLevel == idxContrastHigh)
                    imageContrastLevel = imageContrastLevel + 1;
                    
                    % Play the sound.
                    MakeBeepSound('preset','correct');
                else
                    % Print out it reached the limit.                
                    fprintf('\t You reached the highest contrast! \n');
                    
                    % Play the sound for feedback.
                    MakeBeepSound('preset','incorrect');
                end
            end
            
        elseif strcmp(buttonPress,'left')
            % Show the same contrast level again for next display.
            imageContrastLevel = imageContrastLevel;
            
            % Play the feedback sound.
            numPlaySound = 2;
            for pp = 1:numPlaySound
                MakeBeepSound('preset','correct');
                
            end
        end
        end
    
    % Play sound as feedback when the contrast level was decided.
    numPlaySound = 3;
    for pp = 1:numPlaySound
        MakeBeepSound('preset','correct');
    end
    
    % Print out the contrast level we found.
    thresholdFoundRawLinear(cc) = predefinedContrasts(imageContrastLevel);
    fprintf('Contrast was found at (%.3f) \n', thresholdFoundRawLinear(cc));
    fprintf('Initial contast sensitivity measure is in progress - (%d/%d) \n', cc, nInitialContrasts);
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

% If there is not desired number of test points, add more to match the
% number.
if ~(length(estDomainValidation)==options.nContrastPoints)
    nPredefinedContrasts = length(predefinedContrastsLog);
    nContrastsNeeded = options.nContrastPoints - length(estDomainValidation);
    
    % Get index of predefined contrasts.
    for vv = 1:length(estDomainValidation)
        idxPredefinedContrastsLog(vv) = find(predefinedContrastsLog==estDomainValidation(vv));
    end
    
    idxContrastsNotSelected = setdiff([1:1:nPredefinedContrasts],idxPredefinedContrastsLog);
    
    if any(idxPredefinedContrastsLog,length(estDomainValidation))
        % When the contrast range was set too high.
        idxContrastsNotSelected = sort(idxContrastsNotSelected,'descend');
    else
        % When the contrast range was set too low.
        idxContrastsNotSelected = sort(idxContrastsNotSelected,'aescend');
    end
    
    % Add new contrast points to the pre-selected result.
    idxContrastsToAdd = idxContrastsNotSelected(1:nContrastsNeeded);
    idxContrastsNew = [idxPredefinedContrastsLog idxContrastsToAdd];
    
    % Set the new contrast range.
    estDomainValidation = predefinedContrastsLog(idxContrastsNew);
    estDomainValidation = unique(estDomainValidation);
else
    nContrastsNeeded = 0;
end

% Print out the results.
nContrastPointsSelected = length(estDomainValidation);
if (options.verbose)
    fprintf('The number of contrast points selected = (%d/%d) \n', nContrastPointsSelected, options.nContrastPoints);
    fprintf('\t The number of points (%d/%d) over the predefined range, consider using high contrast image set for 18 cpd \n,',...
        nContrastsNeeded, options.nContrastPoints);
end

% Save the raw data in struct.
preExpDataStruct.rawData.testContrast = testContrastCollect;
preExpDataStruct.rawData.rngValues = rngValuesCollect;
preExpDataStruct.rawData.whichDirectionToDisplay = whichDirectionCollect;

preExpDataStruct.thresholdFoundRawLinear = thresholdFoundRawLinear;
preExpDataStruct.thresholdEstLinear = thresholdEstLinear;
preExpDataStruct.estDomainValidationNominalLinear = estDomainValidationNominalLinear;

end
