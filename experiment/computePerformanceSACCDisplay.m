function [correct, flipTime] = computePerformanceSACCDisplay(nullRGBImage,testRGBImage,...
    theSceneTemporalSupportSeconds,theCrossbarTemporalSupportSeconds,...
    testContrast,window,windowRect,options)
% Run one trial of a psychophysical experiment.
%
% Syntax:
%    [correct, flipTime] = computePerformanceSACCDisplay(nullRGBImage,testRGBImage,...
%    theSceneTemporalSupportSeconds,theCrossbarTemporalSupportSeconds,...
%    testContrast,window,windowRect)
%
% Description:
%     Run one trial of a psychophysical experiment and return correct or
%     incorrect.  The trial is TAFC.  The two stimuli are nullRGBImage and
%     testRGBImage.  Subject is correct if he/she chooses testRGBImage.
%
%     It can be executed both with and without using PTB. In the
%     experiment, we will use the PTB version, but we can also run this by
%     displaying the images on the figure if the option 'simulation' set to
%     true.
%
% Inputs:
%     nullRGBImage                      - Null RGB image as a reference for
%                                         evaluation, which means no contrast
%                                         image in SACC project.
%     testRGBImage                      - Test RGB image to be compared with
%                                         the null image.
%     theSceneTemporalSupportSeconds    - Temporal support vector for the
%                                         test stimuli scene sequences in
%                                         second unit.
%     theCrossbarTemporalSupportSeconds - Temporal support vector for the
%                                         cross fixation image in second
%                                         unit.
%     testContrast                      - Contrast value of the input image.
%                                         This is only used when getting an
%                                         automatic response based on the
%                                         probability estimated from the
%                                         psychometric function with certain
%                                         parameters.
%    window                             - PTB window for opened screen.
%    windowRect                         - Rect corresonding to window.
%
% Outputs:
%     correct                           - 1 if correct and 0 if incorrect.
%     flipTime                          - Array of system flip time of
%                                         displaying images.
%
% Optional key/value pairs:
%     runningMode                       - If you are not acutally running
%                                         the experiment, set this to true and
%                                         null and test images will be showed
%                                         on the figure instead of displaying
%                                         it on PTB.
%    expKeyType                         - Default to 'keyboard'. Set either
%                                         'keyboard' or 'gamepad' to get a
%                                         response when running the
%                                         experiment.
%    imageMagnificationFactor           - When simulation was set to true,
%                                         this decides the size of the null
%                                         and test images. If it is set to 1,
%                                         it displays the original image size.
%    beepSound                          - If it is set to true, make a beep
%                                         sound every when an image is
%                                         displaying as an audible cue for
%                                         patients.
%    autoResponse                       - If it is set to true, make an
%                                         automatic response based on the
%                                         probability
%    debugMode                          - Default to false. If it is set to
%                                         true, a gabor image will be
%                                         displayed until subject press a
%                                         button. It is useful to check the
%                                         stimuli when debugging.
%    movieStimuli                       - Default to false. If it is set to
%                                         true, the stimuli will be
%                                         displayed as gradually ramping of
%                                         and off.
%    movieImageDelaySec                 - Default to 0.5. When using the
%                                         option of 'movieStimuli', this
%                                         decides the duration of the
%                                         gradation of the stimuli changes
%                                         in sec.
%    preStimuliDelaySec                 - Default to 0. Make a time delay
%                                         on null stimulus without a
%                                         crossbar before showing the test
%                                         contrast image.
%    verbose                            - Boolean. Default true. Controls
%                                         printout.
%
% See also
%   t_thresholdEngine, t_spatialCsf, computeThresholdTAFC,
%   SACC_runExperiment

% History:
%    10/23/20  dhb                      - Comments.
%    02/02/22  smo                      - Modified it to use in SACC project.
%    02/04/22  smo                      - Removed the input variable
%                                         displayControlStruct which is not
%                                         needed in this version. We can bring
%                                         it back if we want to.
%    02/07/22  smo                      - Added an option getting automatic
%                                         response. Also, made it shorter and
%                                         clearer in overall.
%    03/17/22  smo                      - Now we set the temporal scene
%                                         support time differently for
%                                         stimuli and cross fixation point.
%    04/27/22  smo                      - Added a running option in PTB to
%                                         display the test image either
%                                         vertical or horizontal.
%    07/11/22  smo                      - Added debug mode option to
%                                         keep displaying a stimulus until
%                                         button pressed.
%    07/13/22  smo                      - Added an option to make stimuli
%                                         presentation gradually ramping
%                                         on and off.
%    07/18/22  smo                      - Added an option to make a time
%                                         delay on null image before
%                                         showing the test contrast image.
%    08/04/22  smo                      - Collecting flip time as an
%                                         output.
%    08/05/22  smo                      - Now we make medium images on the
%                                         cosine function for ramping
%                                         on/off stimuli.

%% Set parameters.
arguments
    nullRGBImage
    testRGBImage
    theSceneTemporalSupportSeconds (1,1)
    theCrossbarTemporalSupportSeconds (1,1)
    testContrast
    window (1,1)
    windowRect (1,4)
    options.runningMode = 'simulation'
    options.expKeyType = 'keyboard'
    options.imageMagnificationFactor (1,1) = 1.3
    options.beepSound (1,1) = false
    options.autoResponse = []
    options.debugMode (1,1) = false
    options.movieStimuli (1,1) = false
    options.movieImageDelaySec (1,1) = 0.5
    options.preStimuliDelaySec (1,1) = 0
    options.verbose (1,1) = true
end

%% Displaying the images.
switch (options.runningMode)
    case 'PTB-sequential'
        % It displays the images on the projector using PTB. The null and
        % test images will be presented sequentially.
        %
        % Randomize the displaying order of null and test images.
        displayFirstTest = 1;
        displayFirstNull = 2;
        displayOrders = [displayFirstTest displayFirstNull];
        whichOneToStart = randi(displayOrders);
        
        allImages = {testRGBImage nullRGBImage};
        firstDisplayImage  = allImages{whichOneToStart};
        secondDisplayImage = allImages{setdiff(displayOrders,whichOneToStart)};
        
        % Display the images here.
        %
        % Cross-fixation image before displaying the images. This is useful
        % when audible cue is not working.
        DisplayScreenPattern(window,windowRect,'patternType','crossbar',...
            'patternColor',[0 0 0],'imageBackground',nullRGBImage,'verbose',false);
        
        % Make a time delay.
        WaitSecs(theCrossbarTemporalSupportSeconds);
        
        % First image.
        SetScreenImage(firstDisplayImage, window, windowRect,'verbose',options.verbose);
        
        % Make a beep sound as an audible cue.
        if (options.beepSound)
            MakeBeepSound;
        end
        
        % Make a time delay.
        WaitSecs(theSceneTemporalSupportSeconds);
        
        % Cross-fixation image again.
        DisplayScreenPattern(window,windowRect,'patternType','crossbar',...
            'patternColor',[0 0 0],'imageBackground',nullRGBImage,'verbose',false);
        
        % Make a time delay.
        WaitSecs(theCrossbarTemporalSupportSeconds);
        
        % Second image.
        SetScreenImage(secondDisplayImage, window, windowRect,'verbose',options.verbose);
        
        % Make a beep sound as an audible cue.
        if (options.beepSound)
            MakeBeepSound;
        end
        
        % Make a time delay.
        WaitSecs(theSceneTemporalSupportSeconds);
        
        % Cross-fixation image again.
        DisplayScreenPattern(window,windowRect,'patternType','crossbar',...
            'patternColor',[0 0 0],'imageBackground',nullRGBImage,'verbose',false);
        
        
    case 'PTB-directional'
        % This part will be used for the actual experiment displaying the test
        % image on the projector using PTB. It displays a test image either
        % vertical or horizontal to make an evaluation, so null image is
        % not displayed during the session after the initial presentation.
        
        % Randomize which direction of the image to display.
        displayVertical   = 1;
        displayHorizontal = 2;
        displayDirections = [displayVertical displayHorizontal];
        whichDirectionToDisplay = randi(displayDirections);
        
        % Make a rotation image here.
        switch whichDirectionToDisplay
            case displayVertical
                % If it's vertical image to display, just pass the original
                % image as the images were modulated in vertical patterns.
                displayTestImage = testRGBImage;
            case displayHorizontal
                % If it's horizontal image to display, rotate the original
                % test image 90 degrees.
                imgRotationDeg = 90;
                displayTestImage = imrotate(testRGBImage, imgRotationDeg);
        end
        
        % Display crossbar image.
        %
        % Cross-fixation image before displaying the images.
        flipTimeInitial = DisplayScreenPattern(window,windowRect,'patternType','crossbar',...
            'patternColor',[0 0 0],'imageBackground',nullRGBImage,'verbose',false);
        
        % Display null image without crossbar.
        %
        % This part can control the time delay before presenting the test
        % contrast image. It is optional and default set to zero.
        if (~options.preStimuliDelaySec == 0)
            flipTimeNull = SetScreenImage(nullRGBImage,window,windowRect,...
                'afterFlipTimeDelay',options.preStimuliDelaySec,'verbose',false);
        else
            flipTimeNull = [];
        end
        
        % Display a test image gradually on.
        %
        % We can make stimuli gradually ramping on and off if we want.
        % Strategy used here is to make four more images having different
        % contrasts between null (zero contrast) and test (target contrast)
        % images, and present them sequentially before displaying the test
        % image so that it can be perceived as the target image is appeared
        % gradually.
        %
        % Here we make four medium images between null and test image that
        % has 20, 40, 60, and 80 percent of the contrast of the test image.
        %
        % Now we make the medium images on the cosine function.
        if (options.movieStimuli)
            
            % Set the contrast ratio of the medium images on the cosine
            % function.
            horizontalLocationIntervalRatio = 0.2;
            movieContrastRatio = mat2gray(cos([pi : horizontalLocationIntervalRatio*pi : 2*pi]));
            movieContrastRatio = setdiff(movieContrastRatio, [0 1]);
            nMovieContrastRatio = length(movieContrastRatio);
            
            % Make medium images here.
            for ii = 1:nMovieContrastRatio
                movieMediumImages{ii} = displayTestImage * movieContrastRatio(ii) + nullRGBImage * (1-movieContrastRatio(ii));
            end
            
            % Display medium images here before displaying the test image.
            movieEachImageDelaySec = options.movieImageDelaySec/nMovieContrastRatio;
            for ii = 1:nMovieContrastRatio
                flipTimeMovieOn(ii,1) = SetScreenImage(movieMediumImages{ii},window,windowRect,...
                    'afterFlipTimeDelay',movieEachImageDelaySec,'verbose',false);
            end
        end
        
        % Display test image here.
        flipTimeTest = SetScreenImage(displayTestImage,window,windowRect,...
            'afterFlipTimeDelay',theSceneTemporalSupportSeconds,'verbose',false);
        
        % For debug mode only, keep displaying the stimulus until the button
        % pressed.
        if (options.debugMode)
            switch (options.expKeyType)
                case 'gamepad'
                    numButtonRight = 3;
                    responseDebug = GetGamepadResp2AFC('numButtonB',numButtonRight,'verbose',options.verbose);
            end
        end
        
        % Display stimuli gradually off.
        if (options.movieStimuli)
            for ii = 1:nMovieContrastRatio
                flipTimeMovieOff(ii,1) = SetScreenImage(movieMediumImages{nMovieContrastRatio-ii+1},window,windowRect,...
                    'afterFlipTimeDelay',movieEachImageDelaySec,'verbose',false);
            end
        end
        
        % Display the cross-fixation image again.
        flipTimeFinal = DisplayScreenPattern(window,windowRect,'patternType','crossbar',...
            'patternColor',[0 0 0],'imageBackground',nullRGBImage,'verbose',false);
        
        % Collect all flip time here.
        flipTime = [flipTimeInitial; flipTimeNull; flipTimeMovieOn; flipTimeTest; flipTimeMovieOff; flipTimeFinal];
        
    case 'simulation'
        % This part does not use PTB and just diplay test images side-by-side on
        % the figure and record key stroke responses, which would be helpful to
        % check the sequence of the experiment before running it with the patients.
        % Resize the image as desired here.
        nullRGBImageResize = imresize(nullRGBImage,options.imageMagnificationFactor);
        testRGBImageResize = imresize(testRGBImage,options.imageMagnificationFactor);
        imageSize = size(nullRGBImageResize);
        imageXPixel = imageSize(1);
        imageYPixel = imageSize(2);
        
        % Get the display info here.
        screenSize = get(0,'screensize');
        screenXPixel = screenSize(3);
        screenYPixel = screenSize(4);
        
        % Randomize the location of the images.
        imageSideLeft  = 1;
        imageSideRight = 2;
        imageSides = [imageSideLeft imageSideRight];
        whichSideNullImage = randi([imageSideLeft imageSideRight]);
        whichSideTestImage = setdiff(imageSides,whichSideNullImage);
        
        % Display images here.
        figure; clf;
        
        % Set the position and the size of the images.
        %
        % Note that the (x,y) = (0,0) is the left and bottom end of the screen.
        % We will display the images on the full screen anyways.
        x = 0;
        y = 0;
        width = screenXPixel;
        height = screenYPixel;
        set(figure, 'Position', [x y width height])
        
        % Display null image.
        subplot(1,2,whichSideNullImage); imshow(nullRGBImageResize);
        % Display test image.
        subplot(1,2,whichSideTestImage); imshow(testRGBImageResize);
        
        % Make a beep sound as an audible cue.
        if (options.beepSound)
            MakeBeepSound;
        end
end

%% Getting a key response.
%
% Key press response is saved in a single number based on the ASCII
% allocated number for the keyboards.
% [28 = leftarrow / 29 = rightarrow].
%
% Maybe we want to change this key to press for the case using the PTB
% because it is not very intuitive. Or we can switch to use the
% joy-stick to evaluate.
leftArrow  = 28;
rightArrow = 29;

% Make a loop till patients press either left of right arrow key.
if (isempty(options.autoResponse))
    switch (options.expKeyType)
        % Get a response using a keyboard.
        case 'keyboard'
            metCondition = false;
            nLoopTrial = 1;
            while (nLoopTrial >= 1)
                gettingResponse = waitforbuttonpress;
                response = double(get(gcf,'CurrentCharacter'));
                if (response == leftArrow || response == rightArrow)
                    metCondition = true;
                    break
                else
                    metCondition = false;
                    nLoopTrial = nLoopTrial + 1;
                end
            end
            
            % Get a response using a gamepad.
        case 'gamepad'
            gamepadRespFirst  = 1;
            gamepadRespSecond = 2;
            
            % We will use different gamepad keys depending on how we present
            % the stimuli.
            if (options.debugMode)
                responseGamePad = responseDebug;
            else
                switch (options.runningMode)
                    case 'PTB-sequential'
                        responseGamePad = GetGamepadResp2AFC('verbose',options.verbose);
                    case 'PTB-directional'
                        numButtonRight = 3;
                        responseGamePad = GetGamepadResp2AFC('numButtonB',numButtonRight,'verbose',options.verbose);
                end
            end
            
            % Converts the response here.
            if (responseGamePad == gamepadRespFirst)
                response = leftArrow;
            elseif (responseGamePad == gamepadRespSecond)
                response = rightArrow;
            end
    end
    
    % Convert the key response into a single number either 0 or 1.
    correctResponse   = 1;
    incorrectResponse = 0;
    
    % Match the variable names prior to response conversion.
    if (strcmp(options.runningMode,'PTB-sequential'))
        whichSideTestImage = whichOneToStart;
        imageSideLeft = displayFirstTest;
        imageSideRight = displayFirstNull;
    elseif (strcmp(options.runningMode,'PTB-directional'))
        whichSideTestImage = whichDirectionToDisplay;
        imageSideLeft  = displayVertical;
        imageSideRight = displayHorizontal;
    end
    
    % Convert the numbers here.
    switch whichSideTestImage
        % When test image was shown on the left side (simulation) or first order (PTB).
        case (imageSideLeft)
            response(response == rightArrow) = incorrectResponse;
            response(response == leftArrow)  = correctResponse;
            % When test image was shown on the right side (simulation) or second order (PTB).
        case (imageSideRight)
            response(response == rightArrow) = correctResponse;
            response(response == leftArrow)  = incorrectResponse;
        otherwise
            error('Check if your images are placed in right positions!');
    end
    
    % Check if the response value is either 0 or 1.
    if (~any(response == [correctResponse, incorrectResponse]))
        error('Response value should be either 0 or 1!');
    end
    
else
    % Getting an auto response based on the test contrast level.
    nTrials = 1;
    [responseVec] = options.autoResponse.psiFunc(log10(testContrast),options.autoResponse.psiParams);
    response = binornd(nTrials,responseVec(2));
    if (options.verbose)
        fprintf('Test contrast %0.4f, pCorrect = %0.3f, response %d\n',testContrast,responseVec(2),response);
    end
end

%% Close.
switch (options.runningMode)
    case 'PTB-sequential'
    case 'PTB-directional'
    case 'simulation'
        close all;
end

% Print out the results.
correct = response;

if (~isempty(options.autoResponse))
    correctResponse   = 1;
    incorrectResponse = 0;
end

% Play the sound based on the evaluations.
if (options.beepSound)
    switch correct
        case correctResponse
            MakeBeepSound('preset','correct');
        case incorrectResponse
            MakeBeepSound('preset','incorrect');
    end
end

end
