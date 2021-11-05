function [settings] = IntegersToSettings(integers,options)
% Quantize the passed [0,1] settings values to integer values
%
% Syntax: [settings] = IntegersToSettings(integers,options)
%
% Description:
%    Take integers and convert to settings quantized on [0,1] as real numbers,
%    appropriate number of levels for the display being controlled.
%
% Inputs:
%    integers -                   The integer settings.  Can be any shape.
%                                 Quantization is applied entrywise.
%
% Outputs:
%    settings -                   A matrix of settings values between 0 and
%                                 1. Matrix as same shape as the input
%                                 integers.
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
%    11/05/21  dhb                Wrote it

%% Set parameters.
arguments
    integers
    options.nInputLevels (1,1) = 253
    options.verbose (1,1) = true
end

%% Do the conversion
settings = integers/(options.nInputLevels-1);

end