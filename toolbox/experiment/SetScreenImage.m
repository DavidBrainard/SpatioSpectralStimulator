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
%    01/06/21  smo                Started on it

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
% display images using PTB than the function screen('PutImage').
imageTexture = Screen('MakeTexture', window, image);
Screen('DrawTexture', window, imageTexture, [], windowRect);
Screen('Flip', window);

% Show the verbose message.
if (options.verbose)
    fprintf('Image is now being displayed on the screen...\n');
end

end
