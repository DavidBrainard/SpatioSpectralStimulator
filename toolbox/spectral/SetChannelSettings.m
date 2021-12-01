function [] = SetChannelSettings(theSettings,options)
% Set the projector subprimary settings for SACC project.
%
% Syntax: [] = SetChannelSettings(theSettings)
%
% Description:
%    This sets the projector channel settings for the projector.
%
% Inputs:
%    theSettings -                Channel settings we wish to set.  This
%                                 is a matrix with nPrimaries columns (one for
%                                 each screen primary) and nChannels
%                                 rows, one for each channel primary. These
%                                 values are specified as real numbers
%                                 between 0 and 1, and are assumed to be
%                                 after gamma correction.
%
% Optional key/value pairs:
%    'logicalToPhysical' -        Vector containing logical to physical
%                                 channel numbering for this projector.
%                                 This should have the same length as the
%                                 number of channels for the given
%                                 projector.
%    'nInputLevels' -             Number of channel input levels.
%                                 Default is 253.
%    'nPrimaries' -               Number of projector primaries. Typically 3
%                                 but you never know.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    10/28/21  smo                Started on it
%    11/12/21  dhb,smo            Deleted the feature controlling the
%                                 projector mode. It is now controlled
%                                 when opening the projector screen
%                                 (OpenPlainScreen), and a
%                                 separate function for it will be added.
%    12/01/21  smo                Now logicalToPhysical becomes [0:15] as
%                                 our new projector has 16 independent
%                                 working channels.

%% Set parameters.
arguments
    theSettings
    options.logicalToPhysical = [0:15]
    options.nInputLevels (1,1) = 253
    options.nPrimaries (1,1) = 3
    options.verbose (1,1) = true
end

% Check consistency of passed settings
[m,n] = size(theSettings);
if (m ~= length(options.logicalToPhysical))
    error('Number of subprimary settings passed does not match logicalToPhysical mapping');
end
if (n ~= options.nPrimaries)
    error('Number of columns in settings does not match number of projector primaries');
end

% Convert subprimary settings to integers
channelIntegers = SettingsToIntegers(theSettings,'nInputLevels',options.nInputLevels);

% Set projector current levels as the above settings.
nChannels = size(theSettings,1);
for pp = 1:options.nPrimaries
    for ss = 1:nChannels
        Datapixx('SetPropixxHSLedCurrent', pp-1, options.logicalToPhysical(ss), channelIntegers(ss,pp));
    end
end

end