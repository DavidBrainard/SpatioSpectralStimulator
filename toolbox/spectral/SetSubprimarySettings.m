function [] = SetSubprimarySettings(subprimarySettings,options)
% Set the projector subprimary settings for SACC project.
%
% Syntax: [] = SetSubprimarySettings(subprimarySettings)
%
% Description:
%    This sets up the projector subprimary settings for the projector. This
%    should be used wherever the projector displays an image including
%    calibration, displaying gabor patch, etc.
%
% Inputs:
%    subprimarySettings -         Subprimary settings we wish to set.  This
%                                 is a matrix with nPrimaries columns (one for
%                                 each projector primary) and nSubprimaries
%                                 rows, one for each subprimary. These
%                                 values are specified as real numbers
%                                 between 0 and 1, and are assumed to be
%                                 after gamma correction.
%
% Optional key/value pairs:
%    'logicalToPhysical' -        Vector containing logical to physical
%                                 subprimary numbering for this projector.
%                                 This should have the same length as the
%                                 number of subprimaries for the given
%                                 projector.
%    'nInputLevels' -             Number of subprimary input levels.
%                                 Default is 253.
%    'nPrimaries' -               Number of projector primaries. Typicall 3
%                                 but you never know.
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either to be 'Normal' (true) or
%                                 'Steady-on' (false).
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    10/28/21  smo                Started on it

%% Set parameters.
arguments
    subprimarySettings
    options.logicalToPhysical = [0:7 9:15]
    options.nInputLevels (1,1) = 253
    options.nPrimaries (1,1) = 3
    options.projectorMode (1,1) = true
    options.verbose (1,1) = true
end

%% Set the projector mode.
if (options.projectorMode)
    commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Normal mode (Default)
    unix(commandNormal)
    if (options.verbose)
        disp('Projector is set as Normal mode');
    end
else
    commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit'; % Steady-on mode
    unix(commandSteadyOn)
    if (options.verbose)
        disp('Projector is set as Steady-on mode');
    end
end

% Check consistency of passed settings
[m,n] = size(subprimarySettings);
if (n ~= length(options.logicalToPhysical))
    error('Number of subprimary settings passed does not match logicalToPhysical mapping');
end
if (m ~= options.nPrimaries)
    error('Number of columns in subprimary settings does not match number of projector primaries');
end

% Convert subprimary settings to integers
subprimaryIntegers = SettingsToIntegers(subprimarySettings,'nInputLevels',options.nInputLevels);

% Set projector current levels as the above settings.
nSubprimaries = size(subprimarySettings,2);
for pp = 1:options.nPrimaries
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', pp-1, options.logicalToPhysical(ss), subprimaryIntegers(pp,ss));
    end
end

end