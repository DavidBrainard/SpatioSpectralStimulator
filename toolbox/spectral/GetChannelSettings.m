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
%    verbose                    - Default true. Boolean. Controls plotting
%                                 and printout.

% History:
%    03/18/22  smo              - Started on it

%% Set parameters.
arguments
    options.nPrimaries (1,1) = 3;
    options.verbose (1,1) = true;
end

%% Read the channel settings.
%
% Connect to the projector if it hasn't. It returns non zero if Datapixx has
% been successfully opened for use.
if (~isReady)
    isReady = Datapixx('open');
    isReady = Datapixx('IsReady');
end

% Get channel settings here.
for pp = 1:options.nPrimaries
    channelSettings(:,pp) = Datapixx('GetPropixxHSLedCurrents',pp-1)
end
nChannels = size(channelSettings,1);

% Print out.
if (options.verbose)
    for cc = 1:nChannels
        fprintf('Ch%d: %d    %d    %d \n', cc, ...
            channelSetting(cc,1),channelSetting(cc,2),channelSetting(cc,3));
    end
end
end
