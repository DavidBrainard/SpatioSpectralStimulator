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
%    'pitch' -                    Set the pitch of the beep sound. The
%                                 higher the value is, the higher the beep
%                                 sound becomes.
%    'duration' -                 Set duration of the beep sound in second.
%    'verbose' -                  Boolean. Default true. Controls plotting
%                                 and printout.

% History:
%    01/10/21  smo                Started on it

%% Set parameters.
arguments
    options.pitch (1,1) = 1000
    options.duration (1,1) = 0.2
    options.frequency (1,1) = 65000
    options.verbose (1,1) = true
end

%% Make a sound wave.         
%
% Set the period time matrix based on the input parameters. 
periodTime = [0: 1/options.frequency: options.duration];  

% Make a sound wave form here.
soundWave = sin(2 * pi * options.pitch * periodTime);    

%% Play sound here.
%
% The 'Sound' command can be also used, but it cannot play more than 200
% beeps successively, so we used 'Snd' command here.
%
% Check if Snd function works for the device.
sndStatus = Snd('Open');

if (~sndStatus == 0)
    % Playing beep sound using Snd function.
    Snd('Play', soundWave, options.frequency);
    Snd('Wait');
    Snd('Quiet');
    Snd('Close');
else
    % Another function to play beep sound if Snd function doesn't work. It
    % seems Snd function does not work on the Linux machine for SACC
    % project (as of 3/22/22).
    Beeper;
end

end
