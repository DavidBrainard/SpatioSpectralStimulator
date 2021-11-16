function [theIntegers] = SetPlainScreenSettings(theSettings,window,windowRect,options)
% This displays a plain screen on the projector using Psychtoolbox.
%
% Syntax: [theIntegers] = SetPlainScreenSettings(theSettings,window,windowRect)
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
%    theIntegers -                Integer settings matching the desired
%                                 maximum value.
%
% Optional key/value pairs:
%    'projectorToolboxPath' -     Path to the Vpixx control toolbox.  We
%                                 add this to the Matlab path if it isn't
%                                 already there.  This doesn't need to be
%                                 right if the toolbox is already on the
%                                 path.
%    'nInputLevels' -             This number controls the maximum value
%                                 for display input settings. Default to
%                                 256. This would help to better control
%                                 quantization displaying the colors on the
%                                 projector.
%    'verbose' -                  Boolean. Default true. Controls the printout.
%
% See also: OpenPlainScreen, CloseScreen.

% History:
%    10/28/21  smo                Started on it.
%    11/02/21  smo                Add the feature switching the display
%                                 input range from [0,1] to [0,255].
%    11/09/21  smo                Color range is now controlled based on
%                                 the nInputLevels to match up our
%                                 convention.

%% Set parameters.
arguments
    theSettings (3,1) {mustBeInRange(theSettings,0,1,"inclusive")}
    window (1,1)
    windowRect (1,4)
    options.nInputLevels (1,1) = 256
    options.verbose (1,1) = true
end

% Set the PTB display input range from [0,1] to [0,nInputLevels-1]. Default [0,255].
Screen('ColorRange',window,options.nInputLevels-1);

% Scale the settings to match up the working range.
theIntegers = SettingsToIntegers(theSettings,'nInputLevels',options.nInputLevels);

% Check if the integer settings are within the right range.
if any(theIntegers > options.nInputLevels-1)
   error('The integer settings are in the wrong range!');
end

%% Set the color of plain screen on the projector.
Screen('FillRect',window,theIntegers,windowRect);
Screen('Flip', window);
if (options.verbose)
    fprintf('Projector settings [%.0f, %.0f, %.0f] \n',theIntegers(1),theIntegers(2),theIntegers(3));
end

end