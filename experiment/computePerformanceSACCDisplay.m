function [correct, flipTime, rngValues] = computePerformanceSACCDisplay(nullRGBImage,testRGBImage,...
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
%     rngValues                         - Array of rng values that are the
%                                         seed numbers of random noise
%                                         images. We can retrieve the exact
%                                         same random noise image using this
%                                         number.
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
%    rotateImageDeg                     - Default to zero. Rotate the test
%                                         contrast images as much as this
%                                         number is set. Note that it is
%                                         always 90-deg difference between
%                                         two rotations for 2-AFC
%                                         experiment.
%    addNoiseToImage                    - Default to false. If it is set to
%                                         true, add noise to image. We
%                                         added this part to minimize the
%                                         artifacts when we see the image
%                                         on SACCSFA.
%    addFixationPointImage              - Default to false. Add a fixation
%                                         point at the center of the image
%                                         when it sets to true.
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
%    08/08/22  smo                      - Added an option to display a
%                                         small focusing point on the test
%                                         contrast images to minimize the
%                                         artifacts.
%    08/17/22  smo                      - Erased the parts that are not
%                                         used in the main experiment. Also
%                                         now we separate functions for
%                                         making and displaying the image
%                                         texture.
%    08/22/22  smo                      - Added an option adding noise to
%                                         image.
%    08/29/22  smo                      - Added an option to rotate the
%                                         images.
%    09/08/22  smo                      - We print out the seed number for
%                                         random noise images (rngValues).

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
    options.beepSound (1,1) = false
    options.autoResponse = []
    options.debugMode (1,1) = false
    options.rotateImageDeg (1,1) = 0
    options.addNoiseToImage (1,1) = false
    options.addFixationPointImage (1,1) = false
    options.movieStimuli (1,1) = false
    options.movieImageDelaySec (1,1) = 0.5
    options.preStimuliDelaySec (1,1) = 0
    options.verbose (1,1) = true
end

%% Set the test image and make a rotation if needed.
%
% This part will be used for the actual experiment displaying the test
% image on the projector using PTB. It displays a test image either
% vertical or horizontal to make an evaluation, so null image is
% not displayed during the session after the initial presentation.
%
% Randomize which direction of the image to display.
displayVertical   = 1;
displayHorizontal = 2;
displayDirections = [displayVertical displayHorizontal];
whichDirectionToDisplay = randi(displayDirections);

% Make a rotation on image.
rotationImageType = 'crop';

switch whichDirectionToDisplay
    case displayVertical
        % If it's vertical image to display, just pass the original
        % image unless you want to rotate it.
        if (~options.rotateImageDeg == 0)
            displayTestImage = fliplr(imrotate(testRGBImage, +options.rotateImageDeg, rotationImageType));
        else
            displayTestImage = fliplr(testRGBImage);
        end
    case displayHorizontal
        % If it's horizontal image to display, rotate the original
        % test image 90 degrees. However, you can rotate more or less.
        rotateImageHorizontalDeg = 90;
        displayTestImage = fliplr(imrotate(testRGBImage, rotateImageHorizontalDeg+options.rotateImageDeg, rotationImageType));
end

% Fill out the cropped part of the image with null image as background.
%
% When we make a rotation on the image besides multiples of 90-deg, it
% makes the image cropped by leaving the rest part of the image in black.
% So, here we fill out the area with the same pixel as null image so
% that the image will look natural.
pixelCroppedBlack = 0;
pixelNullImage = squeeze(nullRGBImage(1,1,:));
nPrimaries = size(nullRGBImage,3);

for pp = 1:nPrimaries
    displayTestImage(:,:,pp) = changem(displayTestImage(:,:,pp), pixelNullImage(pp), pixelCroppedBlack);
end

%% Make displaying image texture.
%
% Here we make all displaying images in PTB texture so that we can save
% time when displaying each image and also minimize the frame break-up.
%
% Cross-bar image.
[imageTextureCrossbar imageWindowRect rngValCrossbar] = MakeImageTexture(nullRGBImage, window, windowRect, ...
    'addNoiseToImage', options.addNoiseToImage, 'addFixationPoint', true, 'verbose', false);

% Null image without cross-bar.
[imageTextureNull imageWindowRect rngValNull] = MakeImageTexture(nullRGBImage, window, windowRect, ...
    'addNoiseToImage', options.addNoiseToImage, 'addFixationPoint', false, 'verbose', false);

% Ramping on/off medium images.
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
    
    % Set the contrast ratio offlipTime the medium images on the cosine function.
    horizontalLocationIntervalRatio = 0.2;
    movieContrastRatio = mat2gray(cos([pi : horizontalLocationIntervalRatio*pi : 2*pi]));
    movieContrastRatio = setdiff(movieContrastRatio, [0 1]);
    nMovieContrastRatio = length(movieContrastRatio);
    
    % Make medium images here.
    for pp = 1:nMovieContrastRatio
        movieMediumImageTemp = displayTestImage * movieContrastRatio(pp) + nullRGBImage * (1-movieContrastRatio(pp));
        [imageTextureMovie(pp) imageWindowRect rngValMovie(pp,1)] = MakeImageTexture(movieMediumImageTemp, window, windowRect, ...
            'addNoiseToImage', options.addNoiseToImage, 'addFixationPoint', options.addFixationPointImage, 'verbose', false);
    end
end

% Test contrast image.
[imageTextureTest imageWindowRect rngValTest] = MakeImageTexture(displayTestImage, window, windowRect, ...
    'addNoiseToImage', options.addNoiseToImage, 'addFixationPoint', options.addFixationPointImage, 'verbose', false);

% Collect rng values here.
%
% These are seed numbers for random noise image so that we can retrieve the
% noise images that we used.
rngValues = struct('rngValCrossbar',rngValCrossbar,'rngValNull',rngValNull, ...
    'rngValMovie',rngValMovie,'rngValTest',rngValTest);

%% Displaying the image texture.
%
% Cross-fixation image before displaying the images.
flipTimeInitial = FlipImageTexture(imageTextureCrossbar, window, imageWindowRect, 'verbose', false);

% Display null image without crossbar if you want.
%
% This part can control the time delay before presenting the test
% contrast image. It is optional and default set to zero.
if (~options.preStimuliDelaySec == 0)
    flipTimeNull = FlipImageTexture(imageTextureNull, window, imageWindowRect,...
        'afterFlipTimeDelay',options.preStimuliDelaySec,'verbose',false);
else
    flipTimeNull = [];
end

% Display a test image gradually on.
if (options.movieStimuli)
    movieEachImageDelaySec = options.movieImageDelaySec/nMovieContrastRatio;
    for pp = 1:nMovieContrastRatio
        flipTimeMovieOn(pp,1) = FlipImageTexture(imageTextureMovie(pp), window, imageWindowRect,...
            'afterFlipTimeDelay', movieEachImageDelaySec, 'verbose', false);
    end
end

% Display test image here.
flipTimeTest = FlipImageTexture(imageTextureTest, window, imageWindowRect, ...
    'afterFlipTimeDelay', theSceneTemporalSupportSeconds, 'verbose', false);

% For debug mode only, keep displaying the stimulus until the button
% pressed.
if (options.debugMode)
    switch (options.expKeyType)
        case 'gamepad'
            if (options.rotateImageDeg == 0)
                numButtonUp    = 4;
                numButtonRight = 3;
                responseDebug  = GetGamepadResp2AFC('numButtonA', numButtonUp, 'numButtonB',numButtonRight,'verbose',options.verbose);
            else
                numButtonLeft  = 1;
                numButtonRight = 3;
                responseDebug  = GetGamepadResp2AFC('numButtonA', numButtonLeft, 'numButtonB',numButtonRight,'verbose',options.verbose);
            end
    end
end

% Display stimuli gradually off.
if (options.movieStimuli)
    for pp = 1:nMovieContrastRatio
        flipTimeMovieOff(pp,1) = FlipImageTexture(imageTextureMovie(nMovieContrastRatio-pp+1), window, imageWindowRect,...
            'afterFlipTimeDelay', movieEachImageDelaySec, 'verbose', false);
    end
end

% Display the cross-fixation image again.
flipTimeFinal = FlipImageTexture(imageTextureCrossbar, window, imageWindowRect, ...
    'verbose', true);

% Collect all flip time here.
flipTime = [flipTimeInitial; flipTimeNull; flipTimeMovieOn; flipTimeTest; flipTimeMovieOff; flipTimeFinal];

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
                    case 'PTB-directional'
                        if (options.rotateImageDeg == 0)
                            numButtonUp    = 4;
                            numButtonRight = 3;
                            responseGamePad = GetGamepadResp2AFC('numButtonA', numButtonUp, 'numButtonB',numButtonRight,'verbose',options.verbose);
                        else
                            numButtonLeft = 1;
                            numButtonRight = 3;
                            responseGamePad = GetGamepadResp2AFC('numButtonA', numButtonLeft, 'numButtonB',numButtonRight,'verbose',options.verbose);
                        end
                        
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
    whichSideTestImage = whichDirectionToDisplay;
    imageSideLeft  = displayVertical;
    imageSideRight = displayHorizontal;
    
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
%
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
