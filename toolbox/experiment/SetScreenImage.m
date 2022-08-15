function [flipTime] = SetScreenImage(image,window,windowRect,options)
% This displays images on the screen using PTB.
%
% Syntax:
%    [flipTime] = SetScreenImage(image,window,windowRect)
%
% Description:
%    This is to display test images in the experiment for the SACC project
%    using PTB.
%
% Inputs:
%    image -                      Test images to display on the screen.
%                                 This should be in a image format, not a
%                                 cal format. For example, 512 x 512 x 3
%                                 in double.
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Outputs:
%    flipTime                     System flip time of displaying images.
%
% Optional key/value pairs:
%    timeDelay                    Default to 0. Make a time delay before
%                                 making a flip of the image.
%    addFixationPointImage        Default to false. If it is set to true,
%                                 add a fixation point at the center of the
%                                 image. This is useful when you want to
%                                 add one on the stimuli for SACC project.
%    preFlipTimeDelay             Default to 0. Make a time delay before
%                                 the flip of the screen. Unit in seconds.
%    afterFlipTimeDelay           Default to 0. Make a time delay after the
%                                 flip of the screen. Unit in seconds.
%    verbose                      Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    01/06/22  smo                Started on it
%    02/07/22  smo                Now the image is displayed in its size
%                                 at the center of the screen.
%    07/19/22  dhb, smo           We take flip time as an output.
%    08/04/22  smo                Added an option to make time delay before
%                                 and after the flip of the screen.
%    08/10/22  smo                Added an option to add fixation point at
%                                 the center of the image.

%% Set parameters.
arguments
    image
    window (1,1)
    windowRect (1,4)
    options.addFixationPointImage (1,1) = false
    options.preFlipTimeDelay (1,1) = 0
    options.afterFlipTimeDelay (1,1) = 0
    options.verbose (1,1) = true
end

%% Make an image in PTB texture format.
%
% Convert the image format to uint8.
if (class(image) == 'double')
    image = im2uint8(image);
elseif (class(image) == 'uint8')
    image = image;
else
    error('Input image should be in the format either double or uint8');
end

%% Add fixation point at the center of image if you want.
if (options.addFixationPointImage)
    fixPatternType = 'line';
    fixPatternColor = [0 0 0];
    fixSizePixel = 8;
    fixPatternWidth = 5;
    
    image = AddFixPointImage(image, 'patternType', fixPatternType, 'patternColor',fixPatternColor, ...
        'patternSize', fixSizePixel, 'patternWidth', fixPatternWidth);
end

%% Set the size of the image in PTB texture. 
%
% This is faster and more flexible way to display images using PTB than the
% function screen('PutImage'). We will display the images at the center of
% the screen in its image size.
%
% Make image texture.
imageTexture = Screen('MakeTexture', window, image);

% Make image windowRect for placing it at the center of the screen.
centerScreen = [windowRect(3) windowRect(4)] * 0.5;
imageSizeHalf = [size(image,1) size(image,2)] * 0.5;
imageWindowRect = [centerScreen(1)-imageSizeHalf(1) centerScreen(2)-imageSizeHalf(2) ...
    centerScreen(1)+imageSizeHalf(1) centerScreen(2)+imageSizeHalf(2)];

%% Flip screen and display the image here.
%
% Draw a texture here.
Screen('DrawTexture', window, imageTexture, [], imageWindowRect);

% You can make a time delay before and after the flip of the screen. The
% time is corrected using the function not to break the frame.
%
% Make pre-flip time delay.
if (~isempty(options.preFlipTimeDelay))
    preFlipTimeDelay = MatchScreenFrameTime(options.preFlipTimeDelay);
    WaitSecs(preFlipTimeDelay);
end

% Flip happens here.
flipTime = Screen('Flip', window);

% % Uncomment to test raw speed using flip
% nFrames = 100;
% testTimes = zeros(nFrames,1);
% for ii = 1:nFrames
%     testTimes(ii) = Screen('Flip',window);
% end
% hist(diff(testTimes));

% Show the verbose message if you want.
if (options.verbose)
    fprintf('Image is now being displayed on the screen...\n');
end

% Make after-flip time delay.
if (~isempty(options.afterFlipTimeDelay))
    afterFlipTimeDelay = MatchScreenFrameTime(options.afterFlipTimeDelay);
    WaitSecs(afterFlipTimeDelay);
end

end
