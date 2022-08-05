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
%    It seems this function is only working on Linux when you call Beeper
%    function first. Something in Beeper function opens or connects to the
%    device to make things ready to play. Anyways, this function uses the
%    function 'sound' internally and it works well.
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
%    duration                   - Default to 0.2. Set duration of the beep
%                                 sound in second.
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

%% Set parameters.
arguments
    options.preset = []
    options.frequency (1,1) = 2000
    options.duration (1,1) = 0.05
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
% beepSound = MakeBeep(options.frequency, options.duration);
beepSound = sin(2 * pi * options.frequency * (0:options.duration*options.samplingRate-1)/options.samplingRate);

% Play sound here.
% 
% We changed the function from 'Beeper' to 'sound'.
% Beeper(options.frequency, options.volume, options.duration);
a = audioplayer(beepSound,8192);
% a.StopFcn = @localRemoveAudioObj;
a.play;
WaitSecs(0.25);
clear a


function localRemoveAudioObj(obj, ~)
    % clear the objects audio data
    % to reduce memory usage.
    obj.clearAudioData();
    clear obj
end

end
