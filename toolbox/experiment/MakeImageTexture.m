function [imageTexture, imageWindowRect, rngVal] = MakeImageTexture(image,window,windowRect,options)
% Make a PTB texture of an image.
%
% Syntax:
%    [imageTexture, imageWindowRect] = MakeImageTexture(image,window,windowRect)
%
% Description:
%    This is to make PTB texture of given image. As we want to display
%    test stimuli exactly when we want, we separate the old function
%    (SetScreenImage) into two, building PTB texture (this function) and
%    flipping texture (FlipImageTexture). 
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
%    imageTexture               - PTB texture of the image.
%    imageWindowRect            - Rect corresonding to image texture size.
%    rngVal                     - Seed of the random noise using rng
%                                 fucntion so that we can recreate the
%                                 noise pattern when we want.
%
% Optional key/value pairs:
%    addNoiseToImage            - Default to false. If it is set to true,
%                                 add noise to image. We added this part to
%                                 minimize the artifacts when we see the
%                                 image on SACCSFA.
%    addFixationPointImage      - Default to false. If it is set to true,
%                                 add a fixation point at the center of the
%                                 image. This is useful when you want to
%                                 add one on the stimuli for SACC project.
%    verbose                      Boolean. Default true.  Controls plotting
%                                 and printout.
% See also:
%    SetScreenImage, FlipImageTexture

% History:
%    08/17/22      smo          - Wrote it.
%    08/22/22      smo          - Added an option adding noise to image.
%    09/08/22      smo          - We save the seed number for random noise.

%% Set parameters.
arguments
    image
    window (1,1)
    windowRect (1,4)
    options.addNoiseToImage (1,1) = false
    options.addFixationPointImage = []
    options.verbose (1,1) = true
end

%% Add noise to image if you want.
if (options.addNoiseToImage)
   noiseLevel = 3;
   [image rngVal] = AddNoiseToImage(image,'noiseLevel',noiseLevel);
end

%% Add fixation point at the center of image if you want.
if (~isempty(options.addFixationPointImage))
    % Set the type of the fixation point.
    switch options.addFixationPointImage
        case 'crossbar'
            fixPatternType = 'line';
        case 'circle'
            fixPatternType = 'circle';
    end
    
    fixPatternColor = [0 0 0];
    fixSizePixel = 12;
    fixPatternWidth = 5;
    
    % Add fixation point here.
    image = AddFixPointImage(image, 'patternType', fixPatternType, 'patternColor',fixPatternColor, ...
        'patternSize', fixSizePixel, 'patternWidth', fixPatternWidth);
end

%% Convert image format in uint8.
%
% Convert the image format to uint8.
if (class(image) == 'double')
    image = im2uint8(image);
elseif (class(image) == 'uint8')
    image = image;
else
    error('Input image should be in the format either double or uint8');
end

%% Set the size of the PTB texture. 
%
% Make image texture.
imageTexture = Screen('MakeTexture', window, image);

% Make image windowRect for placing it at the center of the screen.
centerScreen = [windowRect(3) windowRect(4)] * 0.5;
imageSizeHalf = [size(image,1) size(image,2)] * 0.5;
imageWindowRect = [centerScreen(1)-imageSizeHalf(1) centerScreen(2)-imageSizeHalf(2) ...
    centerScreen(1)+imageSizeHalf(1) centerScreen(2)+imageSizeHalf(2)];

%% Show the verbose message if you want.
if (options.verbose)
    fprintf('Image texture has been made successfully!\n');
end

end
