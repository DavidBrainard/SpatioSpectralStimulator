function [correct] = computePerformanceSACCDisplay(nullRGBImage,testRGBImage,...
    theSceneTemporalSupportSeconds,displayControlStruct,options)
% Run one trial of a psychophysical experiment.
%
% Syntax:
%    [correct] = computePerformanceSACCDisplay(nullRGBImage,testRGBIimage,...
%    theSceneTemporalSupportSeconds,displayControlStruct)
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
%     nullRGBImage                   - Null RGB image as a reference for
%                                      evaluation, which means no contrast
%                                      image in SACC project.
%     testRGBImage                   - Test RGB image to be compared with
%                                      the null image.
%     theSceneTemporalSupportSeconds - Temporal support vector (in
%                                      seconds) for scene sequences.
%     displayControlStruct*          - TBD. It is not called in this
%                                      function, maybe we want to delete it
%                                      later on.
%
% Outputs:
%     correct                        - 1 if correct and 0 if incorrect.
%
% Optional key/value pairs:
%     simulation                     - If you are not acutally running
%                                      the experiment, set this to true and
%                                      null and test images will be showed
%                                      on the figure instead of displaying
%                                      it on PTB.
%    imageMagnificationFactor        - When simulation was set to true,
%                                      this decides the size of the null
%                                      and test images. If it is set to 1,
%                                      it displays the original image size.
%    beepSound                       - If it is set to true, make a beep
%                                      sound every when an image is
%                                      displaying as an audible cue for
%                                      patients.
%    verbose                         - Boolean. Default true. Controls
%                                      printout.
%
% See also
%   t_thresholdEngine, t_spatialCsf, computeThresholdTAFC,
%   SACC_runExperiment

% History:
%   10/23/20  dhb                   - Comments.
%   02/02/22  smo                   - Modified it to use in SACC project.

%% Set parameters.
arguments
    nullRGBImage
    testRGBImage
    theSceneTemporalSupportSeconds
    displayControlStruct
    options.simulation (1,1) = true
    options.imageMagnificationFactor (1,1) = 1.5
    options.beepSound (1,1) = false
    options.verbose (1,1) = true
end

%% Running trials - Using PTB.
%
% This part will be used for the actual experiment displaying the test
% image on the projector using PTB.
if (~options.simulation)
    % Open the screen ready.
    initialScreenSettings = [1 1 1];
    [window windowRect] = OpenPlainScreen(initialScreenSettings,'verbose',options.verbose);
    
    % Display the test images here.
    %
    % Null RGB Image.
    SetScreenImage(nullRGBImage, window, windowRect,'verbose',options.verbose);
    % Make a beep sound as an audible cue.
    if (options.beepSound)
        MakeBeepSound;
    end
    
    % Make a time delay before displaying the other image of the
    % pair.
    for dd = 1:theSceneTemporalSupportSeconds;
        pause(1);
    end
    
    % Test RGB Image.
    SetScreenImage(testRGBImage, window, windowRect,'verbose',options.verbose);
    % Make a beep sound as an audible cue.
    if (options.beepSound)
        MakeBeepSound;
    end
    
    if (options.verbose)
        fprintf('Test image is displaying and waiting for the key press... \n');
    end
    
    % Get a key stroke response here.
    gettingResponse = waitforbuttonpress;
    response(tt,ii) = double(get(gcf,'CurrentCharacter'));
    close all;
    if (options.verbose)
        fprintf('     Key input has been received! \n');
    end
end

%% Running trials - Simulation without using PTB.
%
% This part does not use PTB and just diplay test images side-by-side on
% the figure and record key stroke responses, which would be helpful to
% check the sequence of the experiment before running it with the patients.
if (options.simulation)
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
    
    % Display image here.
    figure; clf;
    
    % Set the position and the size of the test image.
    imageFig = figure;
    x = screenXPixel * 0.2;
    y = screenYPixel * 0.2;
    width = imageXPixel * 2;
    height = imageYPixel;
    set(imageFig, 'Position', [x y width height])
    
    % Left side image.
    subplot(1,2,1); imshow(nullRGBImageResize);
    title('Test Image','fontsize',15);
    % Right side image.
    subplot(1,2,2); imshow(testRGBImageResize);
    
    % Make a beep sound as an audible cue.
    if (options.beepSound)
        MakeBeepSound;
    end
    
    if (options.verbose)
        fprintf('Test image is displaying and waiting for the key is pressed... \n');
    end
    
    % Get a response either Yes or No.
    gettingResponse = waitforbuttonpress;
    response = double(get(gcf,'CurrentCharacter'));
    close all;
    if (options.verbose)
        fprintf('     Key input has been received! \n');
    end
end
if(options.verbose)
    fprintf('Test image evalaution complete! \n');
end

%% Convert the response into a single number either 0 or 1.
%
% Key press response is converted based on the ASCII allocated number for
% the keyboards. [28 = leftarrow / 29 = rightarrow].
%
% Current version always displays the test image on the right, so pressing
% the right arrow key receives the correct answer (1). It should be changed
% according to how we display the images later on.
correctResponse   = 1;
incorrectResponse = 0;
leftArrow  = 28;
rightArrow = 29;
response(response == rightArrow) = correctResponse;
response(response == leftArrow)  = incorrectResponse;

correct = response;

%% Close.
if (~options.simulation)
    CloseScreen;
end
