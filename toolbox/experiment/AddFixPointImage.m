function [imageEdit] = AddFixPointImage(image, options)
% Add a fixation point at the center of the image.
%
% Syntax:
%    [imageEdit] = AddFixPointImage(image)
%
% Description:
%    Add a fixation point at the center of the image. You can choose
%    diffrent type of patterns as you like. This is to achieve a better
%    focus at the center of the image during the experiment so that we can
%    minimize the artifacts for SACC project.
%
% Inputs:
%    image                      - Image that you want to add a fixation
%                                 point on.
%
% Outputs:
%    imageEdit                  - Edited image with fixation point
%
% Optional key/value pairs:
%    patternType                - Default to '+'. Choose a desired
%                                 pattern type to display within 
%                                 [+, x, o, *, s (square)].
%    patternColor               - Deafult to black ([0 0 0]). Set color of
%                                 the pattern. Each value should be ranged
%                                 within 0-1.
%    patternSize                - Default to 10. Set the size of the
%                                 fixation point in pixel.
%    verbose                    - Boolean. Default true. Controls
%                                 printout.
%
% See also:
%    N/A

% History:
%   08/08/22 smo                - Wrote it.

%% Set parameters.
arguments
    image
    options.patternType = '+'
    options.patternColor (1,3) = [0 0 0]
    options.patternSize (1,1) = 10
    verbose (1,1) = false
end

%% Set the position of the fixation point.
imageSize = size(image);
imageCenter = imageSize(1) * 0.5;

% Set fixation point on the image here.
imageEdit = insertMarker(image, [imageCenter imageCenter], ...
    options.patternType, 'color',options.patternColor, 'size',options.patternSize);

end