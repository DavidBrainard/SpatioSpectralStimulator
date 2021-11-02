function [] = SetProjectorPlainScreenSettings(theSettings,window,windowRect,options)
% This displays a plain screen on the projector using Psychtoolbox.
%
% Syntax: [] = SetProjectorPlainScreenSettings(theSettings)
%
% Description:
%    This updates the settings for the projector when it is in plain screen
%    mode.
%
% Inputs:
%    theSettings -                Desired color to display on the
%                                 projector. This should be in format of
%                                 3x1 matrix, and each entry should be in
%                                 the range from 0 (black) to 1 (white).
%                                 Each row respectively matches 
%                                 red, green, and blue channels.
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    'projectorToolboxPath' -     Path to the Vpixx control toolbox.  We
%                                 add this to the Matlab path if it isn't
%                                 already there.  This doesn't need to be
%                                 right if the toolbox is already on the
%                                 path.
%    'verbose' -                  Boolean. Default true.  Controls the printout.
%
% See also: OpenProjectorPlainScreen, CloseProjectorScreen.

% History:
%    10/28/21  smo                Started on it

%% Set parameters.
arguments
    theSettings (3,1) {mustBeInRange(projectorDisplayColor,0,1,"inclusive")}
    options.verbose (1,1) = true
    options.projectorToolboxPath = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx';
end

% Convert to [0,1] to [0,255] here, using SettingsToIntegers().

%% Set the color of plain screen on the projector.
Screen('FillRect',window,theSettings,windowRect);
Screen('Flip', window);
if (options.verbose)
    fprintf('Projector settings [%.2f, %.2f, %.2f] \n',theSettings(1),theSettings(2),theSettings(3));
end

end