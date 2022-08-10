function [] = MakeBeepSound(options)
% This plays a short tone as an audible cue.
%
% Syntax: [] = MakeBeepSound()
%
% Description:
%    This is for playing a short beep sound to play as an audible cue.
%    Especially, it is written for the SACC project to give an alarm for
%    the patients during the experiment.
%
%    Now we are using the function audioplayer to play the sound. We used
%    to use 'Beeper' or 'sound', but both of them either crashed or showed
%    errors when playing the sounds in a loop.
%
%    We make audioplayer object to play the sound and clear it after the
%    play not to stack up the objects. It showed the errors when the
%    objects were collected up to the number of 55.
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
%    frequency                  - Default to 5000. Set the frequency of the
%                                 beep sound. The higher the value is, the
%                                 higher the beep sound becomes.
%    duration                   - Default to 0.25. Set duration of the beep
%                                 sound in seconds. It also controls the
%                                 time delay after making the sound as this
%                                 function did not work without making a
%                                 bit of time delay after playing the
%                                 sound.
%    verbose                    - Default true. Controls printout.

% History:
%    01/10/22  smo              - Started on it
%    03/23/22  smo              - Made it simpler and now using Beeper
%                                 command. Also, made an option to choose
%                                 playing preset sounds.
%    07/28/22  smo              - We switched the function to play the
%                                 sound from Beeper to sound as it crashed
%                                 while running the experiment. It still
%                                 need to be tested if it works fine.
%    08/05/22  dhb, smo         - Now we use the function audioplayer out
%                                 of the sound function, and it is working
%                                 fine without errors when we run it
%                                 multiple times.

%% Set parameters.
arguments
    options.preset = []
    options.frequency (1,1) = 5000
    options.duration (1,1) = 0.1
    options.samplingRate (1,1) = 44100
    options.verbose (1,1) = true
end

%% Play sound here.
%
% You can play preset sound if you want.
if (~isempty(options.preset))
    switch options.preset
        case 'correct'
            options.frequency = 5000;
        case 'incorrect'
            options.frequency = 1200;
    end    
end
 
% Make a beep sound here.
frequencyDefault = 8192;
beepSound = sin(2 * pi * frequencyDefault * (0:options.samplingRate-1)/options.samplingRate);

% Play sound here.
a = audioplayer(beepSound, options.frequency);
a.play;
WaitSecs(MatchScreenFrameTime(options.duration));
clear a;

end
