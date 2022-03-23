function [] = MakeBeepSound(options)
% This plays a short tone as an audible cue.
%
% Syntax: [] = MakeBeepSound(options)
%
% Description:
%    This is for playing a short beep sound to play as an audible cue.
%    Especially, it is written for the SACC project to give an alarm for
%    the patients during the experiment.
%
% Inputs:
%    N/A
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    preset                     - Default to empty. Play the preset sound.
%                                 It can be set either 'correct' or
%                                 'incorrect' to give an audible cue for
%                                 each evaluation in psychophysical
%                                 experiment.
%    frequency                  - Default to 400. Set the frequency of the
%                                 beep sound. The higher the value is, the
%                                 higher the beep sound becomes.
%    volume                     - Default to 0.4. It can be set within the
%                                 range of 0-1.
%    duration                   - Default to 0.2. Set duration of the beep
%                                 sound in second.
%    verbose                    - Default true. Controls printout.

% History:
%    01/10/22  smo              - Started on it
%    03/23/22  smo              - Made it simpler and now using Beeper
%                                 command. Also, made an option to choose
%                                 playing preset sounds.

%% Set parameters.
arguments
    options.preset = []
    options.frequency (1,1) = 400
    options.volume (1,1) = 0.4
    options.duration (1,1) = 0.2
    options.verbose (1,1) = true
end

%% Play sound here.
%
% You can play preset sound if you want.
if (~isempty(options.preset))
    switch options.preset
        case 'correct'
            options.frequency = 'high';
        case 'incorrect'
            options.frequency = 'low';
    end    
end

% Play sound here.
Beeper(options.frequency, options.volume, options.duration);

end
