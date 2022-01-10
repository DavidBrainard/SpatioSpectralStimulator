function [] = MakeBeepSound(pitch,duration)
%
% This plays a short tone as an audible cue.
%
% Sample freq in Hz (1000-65535 Hz)
% It was 8192 when using SOUND command at the end of this code
fs=65000;  
           
if (nargin == 0)
   pitch = 1000;          
   duration = [0:1/fs:.2];   
elseif (nargin == 1)
   duration = [0:1/fs:.2];   
elseif (nargin == 2)
   duration = [0:1/fs:duration];  
end

%one possible wave form
wave=sin(2*pi*pitch*duration);    

% Play sound.
%
% SOUND(wave,fs);
% however, SOUND command cannot play more than 200 beeps successively
% and an error message will show as below
% %??? Unable to write into sound device
% Error in ==> C:\MATLAB6p5\toolbox\matlab\audio\private\playsnd.dll
% Error in ==> C:\MATLAB6p5\toolbox\matlab\audio\sound.m 
% On line 36  ==> playsnd(y,fs,bits);
Snd('Play',wave,fs); 
Snd('Wait');
Snd('Quiet');
Snd('Close');