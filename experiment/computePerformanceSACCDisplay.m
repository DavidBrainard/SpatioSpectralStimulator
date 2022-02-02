function [correct] = computePerformanceSACCDisplay(nullRGBImage,testRGBIimage,...
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
% Inputs:
%     nullRGBImage             -
%     testRGBImage             -
%     temporalSupport          - Temporal support vector (in seconds) for
%                                scene sequences.
%     displayControlStruct     -
%
% Outputs:
%     correct                  - 1 ifcorrect and 0 if incorrect.
%
% Optional key/value pairs:
%     simulation               - If you are not acutally running the
%                                experiment, set this to true and image
%                                figure will show up on the screen instead
%                                of displaying it on PTB. It will be useful
%                                to check and debug.
%    verbose                   - Boolean. Default true. Controls
%                                printout.
%
% See also
%   t_thresholdEngine, t_spatialCsf, computeThresholdTAFC

% History:
%   10/23/20  dhb             - Comments.
%   02/02/22  smo             - Modifying it to use in SACC project.

%% Set parameters.
arguments
    nullRGBImage
    testRGBIimage
    theSceneTemporalSupportSeconds
    displayControlStruct
    options.simulation (1,1) = true
    options.verbose (1,1) = true
end

%% Running trials using PTB.
%
% This part will be used for the actual experiment displaying the test
% image on the projector using PTB.
if (~options.simulation)
    % Open the screen ready.
    initialScreenSettings = [1 1 1];
    [winodw windowRect] = OpenPlainScreen(initialScreenSettings,'verbose',options.verbose);
    
    % Display the test images here.
    %
    % Null RGB Image.
    DisplayImagePTB(nullRGBImage, window, windowRect);
    MakeBeep
    % Test RGB Image.
    DisplayImagePTB(testRGBIimage, window, windowRect);
    
    if (options.verbose)
        fprintf('Test image %d - trial %d is displaying and waiting for the key is pressed... \n',ii,tt);
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
% Figure and record key stroke responses, which would be helpful to check
% the sequence of the experiment before running it with the patients.
if (options.simulation)
    % Resize the image as desired here.
    imageMagnificationFactor = 1.5;
    image = imresize(image,imageMagnificationFactor);
    imageSize = size(image);
    imageXPixel = imageSize(1);
    imageYPixel = imageSize(2);
    
    % Get the display info here.
    screenSize = get(0,'screensize');
    screenXPixel = screenSize(3);
    screenYPixel = screenSize(4);
    
    % Display test image.
    %
    % Now it just displays an image, but this part will be substituted with a
    % separate function displaying the image using PTB later on.
    
    figure; clf;
    
    % Set the position and the size of the test image.
    imageFig = figure;
    x = screenXPixel*0.2;
    y = screenYPixel*0.2;
    width = imageXPixel*2;
    height = imageYPixel;
    set(imageFig, 'Position', [x y width height])
    
    % Left side image.
    subplot(1,2,1); imshow(image);
    title(append('Test Image ',num2str(ii),' - Trial ',num2str(tt)),'fontsize',15);
    
    % Right side image.
    subplot(1,2,2); imshow(image);
    
    if (options.verbose)
        fprintf('Test image %d - trial %d is displaying and waiting for the key is pressed... \n',ii,tt);
    end
    
    % Get a response either Yes or No.
    %
    % It can be also done by using 'ginput'. But, here we used the function
    % waitforbuttonpresss.
    %
    % Following is the ASCII allocated number for the keyboards.
    %
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    gettingResponse = waitforbuttonpress;
    response(tt,ii) = double(get(gcf,'CurrentCharacter'));
    close all;
    
    if (options.verbose)
        fprintf('     Key input has been received! \n');
    end
end
if(options.verbose)
    fprintf('Test image %d evalaution complete! \n',ii);
end

%% Show the results.
%
% Convert the response into 0 / 1
leftArrow  = 28;
rightArrow = 29;
response(response == leftArrow)  = 0;
response(response == rightArrow) = 1;

correct = response;

if (~options.simulation)
    CloseScreen;
end

end