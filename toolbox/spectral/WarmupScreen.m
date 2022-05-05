function [] = WarmupScreen(options)
% This is to warm up the projector.
%
% Syntax:
%    [] = WarmupScreen()
%
% Description:
%    This warms up the projector. It opens the plain screen on the
%    DMD and set LED channels to turn on.
%
% Inputs:
%    N/A
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    nScreenPrimaries           - Default to 3. The number of primaries in
%                                 the viewing device primaries.
%    nChannels                  - Default to 16. The number of
%                                 independently working channels per each
%                                 screen primary. The default is set based
%                                 on the projector for SACC which has 16
%                                 LED channels.
%    warmupTimeMin              - Default to blank. The waiting time for the
%                                 screen to warm up in minute unit. The 20
%                                 minutes is recommended for the projector
%                                 in SACC project. The reason the default
%                                 is set to blank is that it can sometimes
%                                 be difficult to kill the timer whenever
%                                 you want to stop warming up under the
%                                 timer counts.
%    verbose                    - Default true. Boolean. Controls printout.

% History:
%    03/23/22  smo              - Started on it
%    05/05/22  smo              - Added an option to choose the projector
%                                 mode.

%% Set parameters.
arguments
    options.nScreenPrimaries (1,1) = 3
    options.nChannels (1,1) = 16
    options.projectorMode = true
    options.warmupTimeMin = []
    options.verbose (1,1) = true
end

%% Open plain screen and set projector settings for warming up the projector.
%
% Open screen here.
screenSettings = ones(options.nScreenPrimaries,1);
OpenPlainScreen(screenSettings,'projectorMode',options.projectorMode);

% Set channel settings here.
channelSettings = ones(options.nChannels, options.nScreenPrimaries);
SetChannelSettings(channelSettings);

% Message when things are ready.
if (options.verbose)
    fprintf('Screen is ready to warm up! We will wait for (%d) mintues... \n',options.warmupTimeMin);
end

% Count time for warming up and alarm when it's done.
if (~isempty(options.warmupTimeMin))
    minToSec = 60;
    warmupTimeSec = options.warmupTimeMin * minToSec;
    pause(warmupTimeSec);
    
    if (options.verbose)
        fprintf('Screen has been warmed up for (%d) minutes! \n',options.warmupTimeMin);
    end
end

end
