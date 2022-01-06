function [] = DisplayImagePTB(image,options)
% This displays images on the screen using PTB.
%
% Syntax: [] = DisplayImagePTB(image)
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
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    'screenSettings' -           Initial screen settings to initiate it.
%                                 We eventually display the images, so this
%                                 won't be visually seen unless the image
%                                 fills up the whole window screen.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    01/06/21  smo                Started on it

%% Set parameters.
arguments
    image
    options.screenSettings = [1 1 1]
    options.verbose (1,1) = true
end

%% Open the screen.
[window, windowRect] = OpenPlainScreen(options.screenSettings,'verbose',options.verbose);

%% Display an image on the PTB screen.
imageTexture = Screen('MakeTexture', window, image);
Screen('DrawTexture', window, imageTexture, [], windowRect);
Screen('Flip', window);

end
