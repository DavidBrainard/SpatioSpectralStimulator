function [imageTexture, imageWindowRect] = MakeImageTexture(image,window,windowRect,options)
% Make a PTB texture of an image.
%
% Syntax:
%    [imageTexture] = MakeImageTexture(image,window,windowRect)
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
%    verbose                      Boolean. Default true.  Controls plotting
%                                 and printout.
% See also:
%    SetScreenImage, FlipImageTexture

% History:
%    08/17/22      smo          - Wrote it.

%% Set parameters.
arguments
    image
    window (1,1)
    windowRect (1,4)
    options.addFixationPointImage (1,1) = false
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

%% Make image windowRect for placing it at the center of the screen.
centerScreen = [windowRect(3) windowRect(4)] * 0.5;
imageSizeHalf = [size(image,1) size(image,2)] * 0.5;
imageWindowRect = [centerScreen(1)-imageSizeHalf(1) centerScreen(2)-imageSizeHalf(2) ...
    centerScreen(1)+imageSizeHalf(1) centerScreen(2)+imageSizeHalf(2)];

%% Show the verbose message if you want.
if (options.verbose)
    fprintf('Image texture has been made successfully!\n');
end

end
