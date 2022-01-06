function [] = DisplayImagePTB(image,options)
% DisplayImagePTB
%
% Syntax:
%
% Description:
%
% This code is to display Gabor patch for SACC project. It is useful when
% your computer cannot configure the target display (LED screen in SACC
% project) This code can be used regardless of the configuration state as
% it uses Psychtoolbox
%
% Inputs:
%    integers -                   ddd
%
% Outputs:
%    settings -                   ddd
%
% Optional key/value pairs:
%    'screenSettings' -           Set
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    01/06/21  smo                Started on it

%% Set parameters.
arguments
    image
    options.screenSettings = [1 1 1]
    options.verbose (1,1) = true
end

%% Initialize
OpenPlainScreen(options.screenSettings,'verbose',options.verbose);

%% PTB pre-setup
PsychDefaultSetup(2);

screens = Screen('Screens');
screenNumber = max(screens);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Read the resolution of the display

%% Display a background as plain screen
% Backgroundcolor can be set, but here it is set as the same size between
% the background and the display window, so it will not be shown on the final result
backgroundColor = [255 255 255]; % Range 0-255 (8 bits)
[window, backgroundRect] = Screen('OpenWindow', screenNumber, backgroundColor, windowRect); % Size of the background = Size of the window

%% Display the gabor patch as 'texture'
gaborPatchImage_texture = Screen('MakeTexture', window, image); % Set the gaborPatchImage as 'texture' in PTB
Screen('DrawTexture', window, gaborPatchImage_texture, [], backgroundRect);
Screen('Flip', window);

end
