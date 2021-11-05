function [window, windowRect] = OpenProjectorPlainScreen(initialSettings,options)
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
%    initialSettings -            Desired color to display on the
%                                 projector. This should be in format of
%                                 3x1 matrix, and each entry should be in
%                                 the range from 0 (black) to 1 (white).
%                                 Each row respectively matches
%                                 red, green, and blue channels.
%
% Outputs:
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Optional key/value pairs:
%    'projectorToolboxPath' -     Path to the Vpixx control toolbox.  We
%                                 add this to the Matlab path if it isn't
%                                 already there.  This doesn't need to be
%                                 right if the toolbox is already on the
%                                 path.
%    'screenNum' -                Set which screen to display.
%    'verbose' -                  Boolean. Default true.  Controls the printout.
%
% See also: SetProjectorPlainScreenSettings,
%     MeasureProjectorPlainScreenSettings, CloseProjectorScreen.

% History:
%    10/28/21  smo                Started on it
%    11/03/21  smo                Added the feature of selecting which
%                                 screen to display.

%% Set parameters.
arguments
    initialSettings (3,1) {mustBeInRange(initialSettings,0,1,"inclusive")}
    options.projectorToolboxPath = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx';
    options.screenNum (1,1)
    options.verbose (1,1) = true
end

%% Set the default PTB setting.
PsychDefaultSetup(2);
screens = Screen('Screens');

% Set which screen to display.
if (~isfield(options,'screenNum'))
    screenNumber = max(screens); % Default.
else
    screenNumber = options.screenNum;
end

% Set white as the default background. The screen will be flipped to the
% desired color with the following command.
white = WhiteIndex(screenNumber);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

%% Set the settings
SetProjectorPlainScreenSettings(initialSettings,window,windowRect,'verbose',options.verbose);

%% Also make sure we can talk to the subprimaries.
%
% Make sure Vpixx toolbox is on the path and do our best to add it if not.
thePath = path;
isThere = findstr(thePath,[filesep 'VPixx']);
if (isempty(isThere))
    addpath(genpath(options.projectorToolboxPath));
    isThere = findstr(thePath,[filesep 'VPixx']);
    if (isempty(isThere))
        error('Unable to add VPixx toolbox to path. Figure out why not and fix.');
    end
end

% Connect to the Vpixx projector.
isReady = Datapixx('open');
if (~isReady)
    error('Datapixx call returns error');
end
isReady = Datapixx('IsReady');
if (~isReady)
    error('Datapixx call returns error');
end

end