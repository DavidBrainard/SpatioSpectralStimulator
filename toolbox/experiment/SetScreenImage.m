function [] = SetScreenImage(image,window,windowRect,options)
% This displays images on the screen using PTB.
%
% Syntax: [] = SetScreenImage(image)
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
%    N/A
%
% Optional key/value pairs:
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    01/06/22  smo                Started on it
%    02/07/22  smo                Now the image is displayed in its size 
%                                 at the center of the screen.  

%% Set parameters.
arguments
    image
    window (1,1)
    windowRect (1,4)
    options.verbose (1,1) = true
end

%% Display an image on the PTB screen.
%
% Convert the image format to uint8.
if (class(image) == 'double')
    image = im2uint8(image);
elseif (class(image) == 'uint8')
    image = image;
else
    error('Input image should be in the format either double or uint8');
end

% Display image here as a texture. This is faster and more flexible way to
% display images using PTB than the function screen('PutImage'). We will
% display the images at the center of the screen in its image size.
%
% Make image texture.
imageTexture = Screen('MakeTexture', window, image);

% Make image windowRect for placing it at the center of the screen.
centerScreen = [windowRect(3) windowRect(4)] * 0.5;
imageSizeHalf = [size(image,1) size(image,2)] * 0.5;
imageWindowRect = [centerScreen(1)-imageSizeHalf(1) centerScreen(2)-imageSizeHalf(2) ...
    centerScreen(1)+imageSizeHalf(1) centerScreen(2)+imageSizeHalf(2)];

% Display image here.
Screen('DrawTexture', window, imageTexture, [], imageWindowRect);
Screen('Flip', window);

% Show the verbose message.
if (options.verbose)
    fprintf('Image is now being displayed on the screen...\n');
end

end
