function [] = CloseScreen(options)
% This closes the projector PTB screen.
%
% Syntax: [] = CloseScreen()
%
% Description:
%    This closes the PTB projector screen. This should be used wherever after 
%    the projector is used to display colors, so this is used as a pair with    
%    the function 'OpenPlainScreen'.  
%
% Inputs:
%    N/A
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
    options.verbose (1,1) = true
end

%% Close PTB screen.
sca;
if (options.verbose)
    disp('Projector screen has been closed');
end

end