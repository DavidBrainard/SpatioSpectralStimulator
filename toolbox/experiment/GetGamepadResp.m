function [buttonPressStr] = GetGamepadResp(options)
% Get a button response from gamepad.
%
% Syntax:
%    [buttonPressStr] = GetGamepadResp()
%
% Description:
%    This function gets a response from 4 buttons on the gamepad and returns
%    a response in string based on the selection. 
%
% Inputs:
%    N/A
%
% Output:
%    buttonPress       - The recorded button press corresponding to the
%                        button mapping. It gives the result in string
%                        within 'up', 'down', 'left', or 'right'.
%
% Optional key/value pairs:
%    pause             - Default to 0.5. Make a time delay after button
%                        pressed. Unit in second. 
%    verbose           - Default to true. Print out more status messages.

% History:
%    09/26/22  smo     - Wrote it.

%% Set parameters.
arguments
    options.numButtonUp (1,1) = 4
    options.numButtonDown (1,1) = 1
    options.numButtonLeft (1,1) = 3
    options.numButtonRight (1,1) = 2
    options.pause (1,1) = 0.5
    options.verbose (1,1) = true
end

%% Get a button response here.
%
% Say hello.
if (options.verbose)
    disp('Waiting for a button to be pressed...');
end

% Get the gamepad info.
gamepadIndex = Gamepad('GetNumGamepads');

% Get a button press here till either Up or Down button is pressed.
stateButtonUp = false;
stateButtonDown = false;
stateButtonLeft = false;
stateButtonRight = false;

while (stateButtonUp == false && stateButtonDown == false && stateButtonLeft == false && stateButtonRight == false)
    stateButtonUp    = Gamepad('GetButton', gamepadIndex, options.numButtonUp);
    stateButtonDown  = Gamepad('GetButton', gamepadIndex, options.numButtonDown);
    stateButtonLeft  = Gamepad('GetButton', gamepadIndex, options.numButtonLeft);
    stateButtonRight = Gamepad('GetButton', gamepadIndex, options.numButtonRight);
end

% Convert the logical to number either 1 or 2 based on the button press.
if (stateButtonUp)
    buttonPressStr = 'up';
elseif (stateButtonDown)
    buttonPressStr = 'down';
elseif (stateButtonLeft)
    buttonPressStr = 'left';
elseif (stateButtonRight)
    buttonPressStr = 'right';
else
    error('Button should be pressed within up, down, left, and right key');
end

% Print out which key was pressed.
if (options.verbose)
    fprintf('Gamepad key has been pressed! - (%s)\n', buttonPressStr);
end

% Pause if you want.
WaitSecs(options.pause);

end