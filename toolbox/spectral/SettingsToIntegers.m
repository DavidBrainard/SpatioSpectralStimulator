function [theIntegers] = SettingsToIntegers(theSettings,options)
% Quantize the passed [0,1] settings values to integer values
%
% Syntax: [theIntegers] = SettingsToIntegers(theSettings)
%
% Description:
%    Take 0 to 1 real numbers and quantize to integer input values at the
%    appropriate number of levels for the display being controlled.
%
% Inputs:
%    theSettings -                A matrix of settings values between 0 and
%                                 1.  Can be any shape.  Quantization is
%                                 applied entrywise.
%
% Outputs:
%    theIntegers -                The converted settings.  Matrix as same
%                                 shape as the input settings.
%
% Optional key/value pairs:
%    'nInputLevels' -             Number of subprimary input levels.
%                                 Default is 253. This is the number of
%                                 input levels for our subprimary projector
%                                 system.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    10/29/21  smo                Started on it

%% Set parameters.
arguments
    theSettings
    options.nInputLevels (1,1) = 253
    options.verbose (1,1) = true
end

%% Do the conversion
theIntegers = round(theSettings*(options.nInputLevels-1));

end