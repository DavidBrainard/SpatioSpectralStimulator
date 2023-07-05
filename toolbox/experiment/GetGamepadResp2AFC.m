function [buttonPress] = GetGamepadResp2AFC(options)
% Get a button response from gamepad for 2-AFC task.
%
% Syntax:
%    [buttonPress] = GetGamepadResp2AFC()
%
% Description:
%    This function get a response from 2 buttons on the gamepad and returns
%    a response either 1 or 2 based on the selection. You can set which
%    buttons to use on the gamepad. Default is set to 1 for (Y) button / 2
%    for (A) button.
%
% Inputs:
%    N/A
%
% Output:
%    buttonPress       - The recorded button press corresponding to the button
%                        mapping. It either gives 1 or 2.
%
% Optional key/value pairs:
%    numButtonA        - Default to 4 which is "Y" on Longitech F310. The
%                        button number on the gamepad corresponding to the
%                        first interval.
%    numButtonB        - Default to 2 which is "A" on Longitech F310. The
%                        button number on the gamepad corresponding to the
%                        second interval.
%    pause             - Default to 0. You can set the time delay after the
%                        evaluation. Unit is in second.
%    verbose           - Default to true. Print out more status messages.

% History:
%    01/20/22  MAB     - Started.
%    02/08/22  smo     - Modified it to use for SACC project.
%    04/28/22  smo     - Updated variable names which can be applicable for
%                        the task choosing either vertical or horizontal.

%% Set parameters.
arguments
    options.numButtonA (1,1) = 3
    options.numButtonB (1,1) = 2
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
stateButtonA = false;
stateButtonB = false;
while (stateButtonA == false && stateButtonB == false)
    stateButtonA = Gamepad('GetButton', gamepadIndex, options.numButtonA);
    stateButtonB = Gamepad('GetButton', gamepadIndex, options.numButtonB);
end

% Convert the logical to number either 1 or 2 based on the button press.
pressButtonA = 1;
pressButtonB = 2;
if (stateButtonA == true)
    buttonPress = pressButtonA;
    pressedButton = 'First';
elseif (stateButtonB == true)
    buttonPress = pressButtonB;
    pressedButton = 'Second';
else
    error('Button should be pressed either first(1) or second(2) interval key');
end

% Print out which key was pressed.
if (options.verbose)
    fprintf('(%s) key has been pressed! \n',pressedButton);
end

% Pause if you want.
WaitSecs(options.pause);

end
