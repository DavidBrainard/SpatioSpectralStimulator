function [] = SetScreenMode(options)
% Set the projector either Normal or Steady-on mode.
%
% Syntax: [] = SetScreenMode()
%
% Description:
%    This sets the projector either Normal or Steady-on mode. Projector
%    mode is basically set when opening the projector screen
%    (OpenPlainScreen), but you can use this function in case you
%    want to set the projector mode in the middle of the code.
%
% Inputs:
%    N/A
%
% Optional key/value pairs:
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either to be 'Normal' (true) or
%                                 'Steady-on' (false).
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    11/15/21  smo                Started on it

%% Set parameters.
arguments
    options.projectroMode (1,1) = true
    options.verbose (1,1) = true
end

%% Set the projector mode here.
if (options.projectorMode)
    % Normal mode (Default)
    commandNormal = 'vputil rw 0x1c8 0x0 -q quit';
    unix(commandNormal);
    if (options.verbose)
        disp('Projector is set as Normal mode');
    end
else
    % Steady-on mode
    commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit';
    unix(commandSteadyOn);
    if (options.verbose)
        disp('Projector is set as Steady-on mode');
    end
end

end