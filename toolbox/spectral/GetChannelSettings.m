function [channelSettings] = GetChannelSettings(options)
% Get the current projector channel settings.
%
% Syntax:
%    [channelSettings] = GetChannelSettings()
%
% Description:
%    Get and read out the current projector channel settings. This is
%    just for the purpose of checking which LED channel is turned on per
%    each screen primary.
%
% Inputs:
%    N/A
%
% Optional key/value pairs:
%    nPrimaries                 - Default to 3. Number of screen primaries
%                                 to control.
%    nChannels                  - Default to 16. Number of LED channels per
%                                 each screen primary. The default is based
%                                 on the projector used for SACC project.
%    verbose                    - Default true. Boolean. Controls plotting
%                                 and printout.

% History:
%    03/18/22  smo              - Started on it

%% Set parameters.
arguments
    options.nPrimaries (1,1) = 3;
    options.nChannels (1,1) = 16;
    options.verbose (1,1) = true;
end

%% Read the channel settings.
%
% Connect to the projector if it hasn't. It returns non zero if Datapixx has
% been successfully opened for use.
if (~exist('isReady'))
    isReady = Datapixx('open');
    isReady = Datapixx('IsReady');
end

% Get channel settings here.
%
% Note that there is a command getting currents info for all channels for
% one screen primary ('GetPropixxHSLedCurrents'), but somehow it did not
% work, so we made a loop for all channels.
for pp = 1:options.nPrimaries
    for cc = 1:options.nChannels
        channelSettings(cc,pp) = Datapixx('GetPropixxHSLedCurrent',pp-1,cc-1);
    end
end

% Print out.
if (options.verbose)
    disp('Current channel settings as follows...');
    for cc = 1:options.nChannels
        fprintf('    (Ch%2.d) : %3.f  %3.f  %3.f \n', cc, ...
            channelSettings(cc,1),channelSettings(cc,2),channelSettings(cc,3));
    end
end
end
