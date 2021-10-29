function [] = OpenProjectorPlainScreen(projectorDisplayColor,options)
% This displays a plain screen on the projector using Psychtoolbox.
%
% Syntax: [] = OpenProjectorPlainScreen(setProjectorDisplayColor)
%
% Description:
%    This displays a plain screen with a desired color on the projector
%    screen. This function should be used wherever it needs to display and
%    measure the projector including calibration, additivity check, etc.
%
% Inputs:
%    projectorDisplayColor -      Desired color to display on the
%                                 projector. This should be in format of
%                                 1x3 matrix, and each column should be in
%                                 the range from 0 (black) to 1 (white).
%                                 Each column respectively matches 
%                                 red, green, and blue channels.
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    'verbose' -                  Boolean. Default true.  Controls the printout.
%
% History:
%    10/28/21  smo                Started on it

%% Set parameters.
arguments
    projectorDisplayColor (3,1) {mustBeInRange(projectorDisplayColor,0,1,"inclusive")}
    options.verbose (1,1) = true
end

%% Set the default PTB setting. 
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = max(screens);

% Set white as the default background. The screen will be flipped with the
% desired color with the following command.
white = WhiteIndex(screenNumber);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

%% Set the color of plain screen on the projector.
Screen('FillRect',window,projectorDisplayColor,windowRect);
Screen('Flip', window);
if (options.verbose)
    fprintf('             Projector primary is set to [%d, %d, %d] \n',projectorDisplayColor(1),projectorDisplayColor(2),projectorDisplayColor(3));
end

end