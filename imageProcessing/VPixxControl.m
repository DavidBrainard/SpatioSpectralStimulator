%% VPixx Projector Sample code
% *** Before to start ***
% 1) Make sure the MATLAB is started with sudo command (cf. sudo matlab on terminal)
% 
% 2) Type tbUse('BrainardLabBase') to set the paths to all available toolboxes
%
% 3) Set path to the VPixx toolbox located at the directory: /home/colorlab/Documents/MATLAB/toolboxes/VPixx

%% Initialize
clear; close all; clc;

tbUse('BrainardLabBase')
addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx')); % add VPixx toolbox ('Datapixx') to path (incl. all subfolders)

%% LED Pulse control 
% Run terminal command 'vputil' (case sensitive) / both unix and system work
% Terminal command 'devsel ppx' selects the PROPixx device (Just type this on the command)

% command = 'vputil devsel ppx -q';
% unix(command)

command_primaries = 'vputil rw 0x1c8 0x7 -q quit'; % Set to 'All Primaries on' without typing the command / '-q' is the key and 'quit' command terminates the vputil
unix(command_primaries)

% VPixx Primary settings (type each on the command screen)
% rw 0x1c8 0x1 % (Primary 0 always ON) > Screen goes red 
% rw 0x1c8 0x2 % (Primary 1 always ON) > Screen goes green
% rw 0x1c8 0x4 % (Primary 2 always ON) > Screen goes blue
% rw 0x1c8 0x7 % (All primaries on) > Screen goes brighter with no chromaticity changes (at least visually)
% rw 0x1c8 0x0 % (Set to default) > Back to the first screen

%% Call out the Datapixx function
% addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx')); % add VPixx toolbox ('Datapixx') to path (incl. all subfolders)

isReady = Datapixx('open');
isReady = Datapixx('IsReady'); % Returns non-0 if a Datapixx has been successfully opened for use.

% Set primaryColor and subColor
% 'primaryColor' should be a value of either 0, 1, or 2.
% 'subColor' is value between 0 and 15. / *Note: selecting subColor 8 does not have any effect

%% Read the current info from the projector

isReady = Datapixx('open');
isReady = Datapixx('IsReady');

% Read all primaryColor and subColor channels
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

currents

% current_R = currents(1,:);
% current_G = currents(2,:);
% current_B = currents(3,:);

% Original commands ('currents' does not work sometimes on Linux box)
% current = Datapixx('GetPropixxHSLedCurrent',primaryColor,subColor) % Single 'subColor' channel
% currents = Datapixx('GetPropixxHSLedCurrents',primaryColor) % All 'subColor' channel

%% Control by each primary color channel

isReady = Datapixx('open'); 
isReady = Datapixx('IsReady');

current = 50;

subcolor_R = 13; % 0-15 / 8 is not working
subcolor_G = 6; % 0-15
subcolor_B = 2; % 0-15

Datapixx('SetPropixxHSLedCurrent',0,subcolor_R,current) % R
Datapixx('SetPropixxHSLedCurrent',1,subcolor_G,current) % G
Datapixx('SetPropixxHSLedCurrent',2,subcolor_B,current) % B

% Read the current settings
for j=1:3; % PrimaryColor (0-2)with 
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

currents


%% 1) Control one SubColor 

isReady = Datapixx('open'); 
isReady = Datapixx('IsReady');

subColor = 0; % 0-15
current = 252;

for i=1:3 % subColor
      Datapixx('SetPropixxHSLedCurrent', i-1, subColor, current); 
end

% Read the current settings
for j=1:3; % PrimaryColor (0-2)with 
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

currents


%% 2) Control all subColors within a Primary color

% Control currents
primaryColor = 1; % 0-2
current = 100; % 0-255

for i=1:16 % subColorrw 0x1c8 0x7
      Datapixx('SetPropixxHSLedCurrent', primaryColor, i-1, current); 
end 

% Read the current settings
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

currents

%% 3) Control all primaryColor and subColor

isReady = Datapixx('open');
isReady = Datapixx('IsReady');

% Control currents
current = 0; % 0-252

for j=1:3 % PrimaryColor
   for i=1:16 % subColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current); 
   end
end

% Read the current settings
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

currents


%% 4) Control one SubColor and the others set as 0 

isReady = Datapixx('open');
isReady = Datapixx('IsReady');

current = 0; % 0-252

for j=1:3 % PrimaryColor
   for i=1:16 % subColor
         Datapixx('SetPropixxHSLedCurrent', j-1, i-1, current); 
   end
end

% Read the current settings
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 

% Set the subColor No. for investication
subColor = 5; % 0-15 / 
current = 50;

for i=1:3 % subColor
      Datapixx('SetPropixxHSLedCurrent', i-1, subColor, current); 
end

% Read the current settings
for j=1:3; % PrimaryColor (0-2)filename = append('sub',num2str(num),'.mat'); 
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

currents

%% Original command for setting the currents
% Datapixx('SavePropixxHSLedCurrents'); % Save the present currents
% Datapixx('SetPropixxHSLedCurrent', primaryColor, subColor, current); % Set a single current
% Datapixx('SetPropixxHSLedCurrents', primaryColor, currents); % Set multiple currents


%% Some notes on the Datapixx function
% Datapixx is a MEX file for precise control of the DataPixx device from
%   VPixx Technologies. It has many functions; type "Datapixx" for a list:
%   	Datapixx
%  
%   For explanation of any particular Datapixx function, just add a question
%   mark "?". E.g. for 'Open', try either of these equivalent forms:
%   	Datapixx('Open?')
%   	Datapixx Open?
%  
%   Most of the time you won't use this function directly, but instead use
%   the PsychDataPixx() function vid = videoinput('gentl', 1, 'Mono8');
%   common tasks.
%  
%   For setup of Datapixx video operations, check the online help of
%   PsychImaging(), which has multiple functions for interacting with the
%   Datapixx device.
%  
%   For an overview of demos and other useful helper functions for the DataPixx,
%   type "help DatapixxToolbox".

%% Setup functions:
% isReady = Datapixx('Open');
% isReady = Datapixx('OpenIP' [, ip]);% Call out the Datapixx function
% isReady = Datapixx('open');
% isReady = Datapixx('IsReady'); % Returns non-0 if a Datapixx has been successfully opened for use.
% selectedDevice = Datapixx('SelectDevice' [, deviceType=-1] [, deviceName]);
% isDatapixx = Datapixx('IsDatapixx');
% isDatapixx2 = Datapixx('IsDatapixx2');
% isDatapixx3 = Datapixx('IsDatapixx3');
% isViewpixx = Datapixx('IsViewpixx');rw 0x1c8 0x0
% isViewpixx3D = Datapixx('IsViewpixDatapixx('SavePropixxHSLedCurrents');x3D');
% isPropixxCtrl = Datapixx('IsPropixxCtrl');
% isPropixx = Datapixx('IsPropixx');
% isPropixx = Datapixx('IsTrackpixx');
% Datapixx('Close');
% 
% % Global system information:
% ramSize = Datapixx('GetRamSize');
% firmwareRev = Datapixx('GetFirmwareRev');
% time = Datapixx('GetTime');
% Datapixx('SetMarker');
% marker = Datapixx('GetMarker');
% supplyVoltage = Datapixx('GetSupplyVoltage');
% supplyCurrent = Datapixx('GetSupplyCurrent');
% is5VFault = Datapixx('Is5VFault');
% tempCelcius = Datapixx('GetTempCelcius');
% tempFarenheit = Datapixx('GetTempFarenheit');
% 
% % DAC (Digital to Analog Converter) subsystem:
% dacNumChannels = Datapixx('GetDacNumChannels');
% dacRanges = Datapixx('GetDacRanges');
% Datapixx('SetDacVoltages', channelVoltagePairs);
% dacVoltages = Datapixx('GetDacVoltages');
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteDacBuffer', bufferData [, bufferAddress=0] [, channelList=[0:numChannels-1]]);
% Datapixx('SetDacSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, channelList=0] [, bufferBaseAddress=0] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartDacSchedule');
% Datapixx('StopDacSchedule');
% status = Datapixx('GetDacStatus');
% 
% % ADC (Analog to Digital Converter) subsystem:
% adcNumChannels = Datapixx('GetAdcNumChannels');
% adcRanges = Datapixx('GetAdcRanges');
% adcVoltages = Datapixx('GetAdcVoltages');
% Datapixx('EnableDacAdcLoopback');
% Datapixx('DisableDacAdcLoopback');
% Datapixx('EnableAdcFreeRunning');
% Datapixx('DisableAdcFreeRunning');
% Datapixx('SetAdcSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, channelList=0] [, bufferBaseAddress=4e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartAdcSchedule');
% Datapixx('StopAdcSchedule');
% [bufferData, bufferTimetags, underflow, overflow] = Datapixx('ReadAdcBuffer', numFrames [, bufferAddress]);
% status = Datapixx('GetAdcStatus');
% 
% % DOUT (Digital Output) subsystem:
% doutNumBits = Datapixx('GetDoutNumBits');
% Datapixx('SetDoutValues', bitValues [, bitMask = hex2dec('00FFFFFF')])
% doutValues = Datapixx('GetDoutValues');
% Datapixx('EnableDoutButtonSchedules' [, mode = 0]);
% Datapixx('DisableDoutButtonSchedules');
% Datapixx('EnableDoutBacklightPulse');
% Datapixx('DisableDoutBacklightPulse');
% Datapixx('EnableDoutBlink')
% Datapixx('DisableDoutBlink')
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteDoutBuffer', bufferData [, bufferAddress=8e6]);
% Datapixx('SetDoutSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, bufferBaseAddress=8e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartDoutSchedule');
% Datapixx('StopDoutSchedule');
% Datapixx('EnablePixelMode' [, mode = 0]);
% Datapixx('DisablePixelMode');
% Datapixx('EnableVsyncMode');
% Datapixx('DisableVsyncMode');
% status = Datapixx('GetDoutStatus');
% 
% % DIN (Digital Input) subsystem:
% dinNumBits = Datapixx('GetDinNumBits');
% dinValues = Datapixx('GetDinValues');
% Datapixx('SetDinDataDirection', directionMask);
% Datapixx('SetDinDataOut', dataOut);
% doutValues = Datapixx('GetDinDataOut');
% Datapixx('SetDinDataOutStrength', strength);
% Datapixx('EnableDinDebounce');
% Datapixx('DisableDinDebounce');
% Datapixx('EnableDoutDinLoopback');
% Datapixx('DisableDoutDinLoopback');
% Datapixx('SetDinLog' [, bufferBaseAddress=12e6] [, numBufferFrames=1000]);
% Datapixx('StartDinLog');
% Datapixx('StopDinLog');
% [logData, logTimetags, underflow] = Datapixx('ReadDinLog' [, numFrames]);
% status = Datapixx('GetDinStatus');
% 
% % TOUCHPixx (touch panel) subsystem:
% Datapixx('EnableTouchpixx' [, touchPanelMode=0]);
% Datapixx('DisableTouchpixx');
% coordinates = Datapixx('GetTouchpixxCoordinates');
% Datapixx('SetTouchpixxStabilizeDuration', duration);
% Datapixx('SetTouchpixxLog' [, bufferBaseAddress=12e6] [, numBufferFrames=1000]);
% Datapixx('StartTouchpixxLog');
% Datapixx('StopTouchpixxLog');
% [logCoords, logTimetags, underflow] = Datapixx('ReadTouchpixxLog' [, numFrames]);
% Datapixx('EnableTouchpixxLogContinuousMode');
% Datapixx('DisableTouchpixxLogContinuousMode');
% status = Datapixx('GetTouchpixxStatus');
% 
% % Audio Output subsystem:
% Datapixx('InitAudio');
% Datapixx('SetAudioVolume', volume [, source=0]);
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteAudioBuffer', bufferData [, bufferAddress=16e6]);
% delay = Datapixx('GetAudioGroupDelay', sampleRate);
% Datapixx('SetAudioSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, lrMode=mono] [, bufferBaseAddress=16e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartAudioSchedule');
% Datapixx('StopAudioSchedule');
% status = Datapixx('GetAudioStatus');
% 
% % Microphone input subsystem:
% Datapixx('SetMicrophoneSource', source [, gain]);
% Datapixx('EnableAudioLoopback');
% Datapixx('DisableAudioLoopback');
% delay = Datapixx('GetMicrophoneGroupDelay', sampleRate);
% Datapixx('SetMicrophoneSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, lrMode=mono] [, bufferBaseAddress=20e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartMicrophoneSchedule');
% Datapixx('StopMicrophoneSchedule');
% [bufferData, underflow, overflow] = Datapixx('ReadMicrophoneBuffer', numFrames [, bufferAddress]);
% status = Datapixx('GetMicrophoneStatus');
% 
% % Video subsystem:
% Datapixx('SetVideoMode' [, mode=0]);
% Datapixx('SetVideoGreyscaleMode' [, mode=1]);
% Datapixx('SetVideo% % Setup functions:
% isReady = Datapixx('Open');
% isReady = Datapixx('OpenIP' [, ip]);
% selectedDevice = Datapixx('SelectDevice' [, deviceType=-1] [, deviceName]);
% isReady = Datapixx('IsReady');
% isDatapixx = Datapixx('IsDatapixx');
% isDatapixx2 = Datapixx('IsDatapixx2');
% isDatapixx3 = Datapixx('IsDatapixx3');
% isViewpixx = Datapixx('IsViewpixx');
% isViewpixx3D = Datapixx('IsViewpixx3D');
% isPropixxCtrl = Datapixx('IsPropixxCtrl');
% isPropixx = Datapixx('IsPropixx');
% isPropixx = Datapixx('IsTrackpixx');
% Datapixx('Close');
% 
% % Global system information:
% ramSize = Datapixx('GetRamSize');
% firmwareRev = Datapixx('GetFirmwareRev');
% time = Datapixx('GetTime');
% Datapixx('SetMarker');
% marker = Datapixx('GetMarker');
% supplyVoltage = Datapixx('GetSupplyVoltage');
% supplyCurrent = Datapixx('GetSupplyCurrent');
% is5VFault = Datapixx('Is5VFault');
% tempCelcius = Datapixx('GetTempCelcius');
% tempFarenheit = Datapixx('GetTempFarenheit');
% 
% % DAC (Digital to Analog Converter) subsystem:
% dacNumChannels = Datapixx('GetDacNumChannels');
% dacRanges = Datapixx('GetDacRanges');
% Datapixx('SetDacVoltages', channelVoltagePairs);
% dacVoltages = Datapixx('GetDacVoltages');
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteDacBuffer', bufferData [, bufferAddress=0] [, channelList=[0:numChannels-1]]);
% Datapixx('SetDacSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, channelList=0] [, bufferBaseAddress=0] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartDacSchedule');
% Datapixx('StopDacSchedule');
% status = Datapixx('GetDacStatus');
% 
% % ADC (Analog to Digital Converter) subsystem:
% adcNumChannels = Datapixx('GetAdcNumChannels');
% adcRanges = Datapixx('GetAdcRanges');
% adcVoltages = Datapixx('GetAdcVoltages');
% Datapixx('EnableDacAdcLoopback');
% Datapixx('DisableDacAdcLoopback');
% Datapixx('EnableAdcFreeRunning');
% Datapixx('DisableAdcFreeRunning');
% Datapixx('SetAdcSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, channelList=0] [, bufferBaseAddress=4e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartAdcSchedule');
% Datapixx('StopAdcSchedule');
% [bufferData, bufferTimetags, underflow, overflow] = Datapixx('ReadAdcBuffer', numFrames [, bufferAddress]);
% status = Datapixx('GetAdcStatus');
% 
% % DOUT (Digital Output) subsystem:
% doutNumBits = Datapixx('GetDoutNumBits');
% Datapixx('SetDoutValues', bitValues [, bitMask = hex2dec('00FFFFFF')])
% doutValues = Datapixx('GetDoutValues');
% Datapixx('EnableDoutButtonSchedules' [, mode = 0]);
% Datapixx('DisableDoutButtonSchedules');
% Datapixx('EnableDoutBacklightPulse');
% Datapixx('DisableDoutBacklightPulse');
% Datapixx('EnableDoutBlink')
% Datapixx('DisableDoutBlink')
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteDoutBuffer', bufferData [, bufferAddress=8e6]);
% Datapixx('SetDoutSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, bufferBaseAddress=8e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartDoutSchedule');
% Datapixx('StopDoutSchedule');
% Datapixx('EnablePixelMode' [, mode = 0]);
% Datapixx('DisablePixelMode');
% Datapixx('EnableVsyncMode');
% Datapixx('DisableVsyncMode');
% status = Datapixx('GetDoutStatus');
% 
% % DIN (Digital Input) subsystem:
% dinNumBits = Datapixx('GetDinNumBits');
% dinValues = Datapixx('GetDinValues');
% Datapixx('SetDinDataDirection', directionMask);
% Datapixx('SetDinDataOut', dataOut);
% doutValues = Datapixx('GetDinDataOut');
% Datapixx('SetDinDataOutStrength', strength);
% Datapixx('EnableDinDebounce');
% Datapixx('DisableDinDebounce');
% Datapixx('EnableDoutDinLoopback');
% Datapixx('DisableDoutDinLoopback');
% Datapixx('SetDinLog' [, bufferBaseAddress=12e6] [, numBufferFrames=1000]);
% Datapixx('StartDinLog');
% Datapixx('StopDinLog');
% [logData, logTimetags, underflow] = Datapixx('ReadDinLog' [, numFrames]);
% status = Datapixx('GetDinStatus');
% 
% % TOUCHPixx (touch panel) subsystem:
% Datapixx('EnableTouchpixx' [, touchPanelMode=0]);
% Datapixx('DisableTouchpixx');
% coordinates = Datapixx('GetTouchpixxCoordinates');
% Datapixx('SetTouchpixxStabilizeDuration', duration);
% Datapixx('SetTouchpixxLog' [, bufferBaseAddress=12e6] [, numBufferFrames=1000]);
% Datapixx('StartTouchpixxLog');
% Datapixx('StopTouchpixxLog');
% [logCoords, logTimetags, underflow] = Datapixx('ReadTouchpixxLog' [, numFrames]);
% Datapixx('EnableTouchpixxLogContinuousMode');
% Datapixx('DisableTouchpixxLogContinuousMode');
% status = Datapixx('GetTouchpixxStatus');
% 
% % Audio Output subsystem:
% Datapixx('InitAudio');
% Datapixx('SetAudioVolume', volume [, source=0]);
% [nextBufferAddress, underflow, overflow] = Datapixx('WriteAudioBuffer', bufferData [, bufferAddress=16e6]);
% delay = Datapixx('GetAudioGroupDelay', sampleRate);
% Datapixx('SetAudioSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, lrMode=mono] [, bufferBaseAddress=16e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartAudioSchedule');
% Datapixx('StopAudioSchedule');
% status = Datapixx('GetAudioStatus');
% 
% % Microphone input subsystem:
% Datapixx('SetMicrophoneSource', source [, gain]);
% Datapixx('EnableAudioLoopback');
% Datapixx('DisableAudioLoopback');
% delay = Datapixx('GetMicrophoneGroupDelay', sampleRate);
% Datapixx('SetMicrophoneSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, lrMode=mono] [, bufferBaseAddress=20e6] [, numBufferFrames=maxScheduleFrames]);
% Datapixx('StartMicrophoneSchedule');
% Datapixx('StopMicrophoneSchedule');
% [bufferData, underflow, overflow] = Datapixx('ReadMicrophoneBuffer', numFrames [, bufferAddress]);
% status = Datapixx('GetMicrophoneStatus');
% 
% % Video subsystem:
% Datapixx('SetVideoMode' [, mode=0]);
% Datapixx('SetVideoGreyscaleMode' [, mode=1]);
% Datapixx('SetVideoClut', clut);
% Datapixx('SetVideoClutTransparencyColor', color);
% Datapixx('EnableVideoClutTransparencyColorMode');
% Datapixx('DisableVideoClutTransparencyColorMode');
% Datapixx('SetVideoHorizontalSplit', mode(0=MIRROR,1=SPLIT,2=AUTO));
% Datapixx('SetVideoVerticalStereo', mode(0=NO_STEREO,1=STEREO,2=AUTO));
% Datapixx('SetVideoStereoEye', eye(1=Left,0=Right));
% Datapixx('EnableVideoStereoBlueline');
% Datapixx('DisableVideoStereoBlueline');
% Datapixx('SetVideoStereoVesaWaveform', waveform);
% Datapixx('EnableVideoHorizontalOverlay');
% Datapixx('DisableVideoHorizontalOverlay');
% Datapixx('SetVideoHorizontalOverlayBounds', boundsRect);
% Datapixx('SetVideoHorizontalOverlayAlpha', alphaTable);
% Datapixx('SetVideoPixelSyncLine', rasterLine [, singleLine=1] [, blankLine=1]);
% Datapixx('EnableVideoScanningBacklight');
% Datapixx('DisableVideoScanningBacklight');
% Datapixx('EnableVideoRescanWarning');
% Datapixx('DisableVideoRescanWarning');
% Datapixx('SetVideoBacklightIntensity', intensity);
% Datapixx('EnableVideoLcd3D60Hz');
% Datapixx('DisableVideoLcd3D60Hz');
% pixels = Datapixx('GetVideoLine' [, nPixels=HORIZONTAL_RESOLUTION]);
% status = Datapixx('GetVideoStatus');
% Datapixx('SetVideoConsoleDisplay' [, presetDisposition=0]);
% 
% % PROPixx-specific routines:
% Datapixx('SetPropixxDlpSequenceProgram' [, program=0]);
% Datapixx('EnablePropixxCeilingMount');
% Datapixx('DisablePropixxCeilingMount');
% Datapixx('EnablePropixxRearProjection');
% Datapixx('DisablePropixxRearProjection');
% Datapixx('SetPropixx3DCrosstalk', crosstalk);
% Datapixx('SetPropixx3DCrosstalkLR', crosstalk);
% Datapixx('SetPropixx3DCrosstalkRL', crosstalk);
% Datapixx('EnablePropixxLampLed');
% Datapixx('DisablePropixxLampLed');
% Datapixx('EnableHotspotCorrection');
% Datapixx('DisableHotspotCorrection');
% Datapixx('SetPropixxAwake');
% Datapixx('SetPropixxSleep');
% isAwake = Datapixx('IsPropixxAwake');
% Datapixx('SetPropixxLedMask' [, mask=0]);
% Datapixx('EnablePropixxQuad4x3d');
% Datapixx('DisablePropixxQuad4x3d');
% Datapixx('SetGrayLEDCurrents' [, index=0]);
% Datapixx('SetCustomCalibrationCurrents' [, index=0])
% 
% % PROPixx T-Scope Mode routines
% Datapixx('EnablePropixxTScope');
% Datapixx('DisablePropixxTScope');
% Datapixx('EnablePropixxTScopePrepRequest');
% Datapixx('DisablePropixxTScopePrepRequest');
% Datapixx('WritePropixxTScopePages', pageData, pageIndex [, nPages=1]);
% Datapixx('SetPropixxTScopeSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, startPage=0] [, nPages=maxScheduleFrames]);
% Datapixx('StartPropixxTScopeSchedule');
% Datapixx('StopPropixxTScopeSchedule');
% Datapixx('SetPropixxTScopeMode' [, mode=0]);
% Datapixx('SetPropixxTScopeProgramAddress' [, addr=0]);
% Datapixx('SetPropixxTScopeProgramOffsetPage', offset);
% Datapixx('SetPropixxTScopeProgram', program);
% 
% % PROPixx Gaze Contingent Display (GCD) Mode routines
% Datapixx('EnablePropixxTScopeQuad');
% Datapixx('DisablePropixxTScopeQuad');
% Datapixx('EnablePropixxGcdShift');
% Datapixx('DisablePropixxGcdShift');
% Datapixx('EnablePropixxGcdShiftSubframe');
% Datapixx('DisablePropixxGcdShiftSubframe');
% Datapixx('EnablePropixxGcdShiftHardware');
% Datapixx('DisablePropixxGcdShiftHardware');
% Datapixx('EnableGcdShiftHardwareBridge');
% Datapixx('DisableGcdShiftHardwareBridge');
% Datapixx('SetGcdShiftHardwareMode', mode);
% Datapixx('EnablePropixxSoftwareTestPatternLoad');
% Datapixx('DisablePropixxSoftwareTestPatternLoad');
% Datapixx('SetGcdShiftHardwareTransform', xGain, xOffset, yGain, yOffset);
% Datapixx('SetPropixxSoftwareTestPatternLoadPage', page);
% 
% % PROPixx + DP3 Gaze Contingent Display (GCD) Mode routines
% Datapixx('SetGCDEyePosition', xScreen, yScreen);
% Datapixx('EnableVideoDataToPPx');
% Datapixx('DisableVideoDataToPPx');
% 
% % TRACKPixx (any kind) Functions:
% Datapixx('GetEyeDuringCalibration', xScreen, yScreen);
% [xRawRight yRawRight xRawLeft yRawLeft] = Datapixx('GetEyeDuringCalibrationRaw', xScreen, yScreen);
% Datapixx('FinishCalibration');
% calibrations_coeff = Datapixx('GetCalibrationCoeff');
% [xScreenRight yScreenRight xScreenLeft yScreenLeft xRawRight yRawRight xRawLeft yRawLeft timetag] = Datapixx('GetEyePosition');
% convertedArray = Datapixx('ConvertCoordSysToCartesian', sourceArray, offsetX, scaleX, offsetY, scaleY);
% convertedArray = Datapixx('ConvertCoordSysToCustom', sourceArray, offsetX, scaleX, offsetY, scaleY);
% 
% % TRACKPixx3 only Functions:
% Datapixx('SaveCalibration');
% Datapixx('LoadCalibration');
% Datapixx('ClearCalibration');
% image = Datapixx('GetEyeImage');
% Datapixx('SetLedIntensity', ledIntensity);
% ledIntensity = Datapixx('GetLedIntensity');
% Datapixx('SetExpectedIrisSizeInPixels', IrisSize);
% expectedIrisSize = Datapixx('GetExpectedIrisSizeInPixels');
% pupilSize = Datapixx('GetPupilSizeSimple');
% [ppLeftMajor ppLeftMinor ppRightMajor ppRightMinor] = Datapixx('GetPupilSize');
% [ppLeftX ppLeftY ppRightX ppRightY] = Datapixx('GetPupilCoordinatesInPixels');
% [CRLeftX CRLeftY CRRightX CRRightY] = Datapixx('GetCRCoordinatesInPixels');
% Datapixx('SetupTPxSchedule', [bufferbaseAddress=12e6, numberOfEyeData=60000]);
% Datapixx('StopTPxSchedule');
% Datapixx('StartTPxSchedule');
% [bufferData, underflow, overflow] = Datapixx('ReadTPxData', numFrames);
% status = Datapixx('GetTPxStatus');
% Datapixx('ShowOverlay');
% Datapixx('HideOverlay');
% Datapixx('SetTPxSleep');
% Datapixx('SetTPxAwake');
% Datapixx('EnableSearchLimits');
% Datapixx('DisableSearchLimits');
% Datapixx('ClearSearchLimits');
% Datapixx('SetSearchLimits', leftEye, rightEye);
% [leftEye, rightEye] = Datapixx('GetSearchLimits');
% DEPRECATED: Datapixx('EnableTrackpixxAnalogOutput' ,[eyeNumber=0]);
% Datapixx('EnableTPxAna
% DEPRECATED: Datapixx('DisableTrackpixxAnalogOutput');
% Datapixx('DisableTPxAnalogOut');
% Datapixx('PpSizeCalGetData');
% Datapixx('PpSizeCalGetDataComplete');
% Datapixx('PpSizeCalLinearRegression');
% Datapixx('PpSizeCalClear');
% Datapixx('PpSizeCalSet');
% [slope_r_x, slope_r_y, slope_l_x, slope_l_y] = Datapixx('PpSizeCalGet');
% fov_h = Datapixx('GetHorizontalFOV');
% fov_v = Datapixx('GetVerticalFOV');
% pxSize = Datapixx('GetPixelSize');
% pxDensity = Datapixx('GetPixelDensity');
% Datapixx('SetLens', lens);
% lens = Datapixx('GetLens');
% Datapixx('SetDistance', distance);
% distance = Datapixx('GetDistance');
% [leftFixationFlag, rightFixationFlag] = Datapixx('IsSubjectFixating');
% [leftSaccadeFlag, rightSaccadeFlag] = Datapixx('IsSubjectMakingSaccade');
% Datapixx('SetFixationThresholds' [, maxSpeed=2500] [, minNumberOfConsecutiveSamples=25]);
% Datapixx('SetSaccadeThresholds' [, minSpeed=10000] [, minNumberOfConsecutiveSamples=10]);
% 
% % TRACKPixx /mini Only Functions:
% time = Datapixx('GetTimeTPxMini');
% Datapixx('OpenTPxMini', [screenProportion=80]);
% Datapixx('CloseTPxMini');
% Datapixx('CalibrateTargetTPxMini', targetId);
% Datapixx('FinishCalibrationTPxMini');
% Datapixx('GetEyeImageTPxMini', imageArray);
% targets, targetsCount = Datapixx('InitializeCalibrationTPxMini');
% [bufferData] = Datapixx('ReadTPxMiniData');
% Datapixx('SetScreenProportionTPxMini', proportion]);
% Datapixx('SetupTPxMiniRecording', numberOfFrame);
% status = Datapixx('GetTPxMiniStatus');
% Datapixx('RecordTPxMiniData');
% Datapixx('StopTPxMiniRecording');
% [bufferData] = Datapixx('GetTPxMiniRecordedData');
% [bufferData] = Datapixx('GetTPxMiniLastEyePosition');
% 
% % Reading and writing register cache:
% Datapixx('RegWr');
% Datapixx('RegWrRd');
% Datapixx('RegWrVideoSync');
% Datapixx('RegWrRdVideoSync');
% Datapixx('RegWrPixelSync', pixelSequence [, timeout=255]);
% isTimeout = Datapixx('RegWrRdPixelSync, pixelSequence [, timeout=255]);
% 
% % Miscellaneous Routines:
% Datapixx('StopAllSchedules');
% error = Datapixx('GetError');
% Datapixx('ClearError');
% Datapixx('Reset');
% Datapixx('ResetAll');
% VERSION: 3.8.8370 [13/OCT/2020]
% Uses VPixx Device Server 

% Datapixx('SetVideoClutTransparencyColor', color);
% Datapixx('EnableVideoClutTransparencyColorMode');
% Datapixx('DisableVideoClutTransparencyColorMode');
% Datapixx('SetVideoHorizontalSplit', mode(0=MIRROR,1=SPLIT,2=AUTO));
% Datapixx('SetVideoVerticalStereo', mode(0=NO_STEREO,1=STEREO,2=AUTO));
% Datapixx('SetVideoStereoEye', eye(1=Left,0=Right));
% Datapixx('EnableVideoStereoBlueline');
% Datapixx('DisableVideoStereoBlueline');
% Datapixx('SetVideoStereoVesaWaveform', waveform);
% Datapixx('EnableVideoHorizontalOverlay');
% Datapixx('DisableVideoHorizontalOverlay');
% Datapixx('SetVideoHorizontalOverlayBounds', boundsRect);
% Datapixx('SetVideoHorizontalOverlayAlpha', alphaTable);
% Datapixx('SetVideoPixelSyncLine', rasterLine [, singleLine=1] [, blankLine=1]);
% Datapixx('EnableVideoScanningBacklight');
% Datapixx('DisableVideoScanningBacklight');
% Datapixx('EnableVideoRescanWarning');
% Datapixx('DisableVideoRescanWarning');
% Datapixx('SetVideoBacklightIntensity', intensity);
% Datapixx('EnableVideoLcd3D60Hz');
% Datapixx('DisableVideoLcd3D60Hz');
% pixels = Datapixx('GetVideoLine' [, nPixels=HORIZONTAL_RESOLUTION]);
% status = Datapixx('GetVideoStatus');
% Datapixx('SetVideoConsoleDisplay' [, presetDisposition=0]);
% 
% % PROPixx-specific routines:
% Datapixx('SetPropixxDlpSequenceProgram' [, program=0]);
% Datapixx('EnablePropixxCeilingMount');
% Datapixx('DisablePropixxCeilingMount');
% Datapixx('EnablePropixxRearProjection');
% Datapixx('DisablePropixxRearProjection');
% Datapixx('SetPropixx3DCrosstalk', crosstalk);
% Datapixx('SetPropixx3DCrosstalkLR', crosstalk);
% Datapixx('SetPropixx3DCrosstalkRL', crosstalk);
% Datapixx('EnablePropixxLampLed');
% Datapixx('DisablePropixxLampLed');
% Datapixx('EnableHotspotCorrection');
% Datapixx('DisableHotspotCorrection');
% Datapixx('SetPropixxAwake');
% Datapixx('SetPropixxSleep');
% isAwake = Datapixx('IsPropixxAwake');
% Datapixx('SetPropixxLedMask' [, mask=0]);
% Datapixx('EnablePropixxQuad4x3d');
% Datapixx('DisablePropixxQuad4x3d');
% Datapixx('SetGrayLEDCurrents' [, index=0]);
% Datapixx('SetCustomCalibrationCurrents' [, index=0])
% 
% % PROPixx T-Scope Mode routines
% Datapixx('EnablePropixxTScope');
% Datapixx('DisablePropixxTScope');
% Datapixx('EnablePropixxTScopePrepRequest');
% Datapixx('DisablePropixxTScopePrepRequest');
% Datapixx('WritePropixxTScopePages', pageData, pageIndex [, nPages=1]);
% Datapixx('SetPropixxTScopeSchedule', scheduleOnset, scheduleRate, maxScheduleFrames [, startPage=0] [, nPages=maxScheduleFrames]);
% Datapixx('StartPropixxTScopeSchedule');
% Datapixx('StopPropixxTScopeSchedule');
% Datapixx('SetPropixxTScopeMode' [, mode=0]);
% Datapixx('SetPropixxTScopeProgramAddress' [, addr=0]);
% Datapixx('SetPropixxTScopeProgramOffsetPage', offset);
% Datapixx('SetPropixxTScopeProgram', program);
% 
% % PROPixx Gaze Contingent Display (GCD) Mode routines
% Datapixx('EnablePropixxTScopeQuad');
% Datapixx('DisablePropixxTScopeQuad');
% Datapixx('EnablePropixxGcdShift');
% Datapixx('DisablePropixxGcdShift');
% Datapixx('EnablePropixxGcdShiftSubframe');
% Datapixx('DisablePropixxGcdShiftSubframe');
% Datapixx('EnablePropixxGcdShiftHardware');
% Datapixx('DisablePropixxGcdShiftHardware');
% Datapixx('EnableGcdShiftHardwareBridge');
% Datapixx('DisableGcdShiftHardwareBridge');
% Datapixx('SetGcdShiftHardwareMode', mode);
% Datapixx('EnablePropixxSoftwareTestPatternLoad');
% Datapixx('DisablePropixxSoftwareTestPatternLoad');
% Datapixx('SetGcdShiftHardwareTransform', xGain, xOffset, yGain, yOffset);
% Datapixx('SetPropixxSoftwareTestPatternLoadPage', page);
% 
% % PROPixx + DP3 Gaze Contingent Display (GCD) Mode routines
% Datapixx('SetGCDEyePosition', xScreen, yScreen);
% Datapixx('EnableVideoDataToPPx');
% Datapixx('DisableVideoDataToPPx');
% 
% % TRACKPixx (any kind) Functions:
% Datapixx('GetEyeDuringCalibration', xScreen, yScreClut', clut);en);
% [xRawRight yRawRight xRawLeft yRawLeft] = Datapixx('GetEyeDuringCalibrationRaw', xScreen, yScreen);
% Datapixx('FinishCalibration');
% calibrations_coeff = Datapixx('GetCalibrationCoeff');
% [xScreenRight yScreenRight xScreenLeft yScreenLeft xRawRight yRawRight xRawLeft yRawLeft timetag] = Datapixx('GetEyePosition');
% convertedArray = Datapixx('ConvertCoordSysToCartesian', sourceArray, offsetX, scaleX, offsetY, scaleY);
% convertedArray = Datapixx('ConvertCoordSysToCustom', sourceArray, offsetX, scaleX, offsetY, scaleY);
% 
% % TRACKPixx3 only Functions:
% Datapixx('SaveCalibration');
% Datapixx('LoadCalibration');
% Datapixx('ClearCalibration');
% image = Datapixx('GetEyeImage');
% Datapixx('SetLedIntensity', ledIntensity);
% ledIntensity = Datapixx('GetLedIntensity');
% Datapixx('SetExpectedIrisSizeInPixels', IrisSize);
% expectedIrisSize = Datapixx('GetExpectedIrisSizeInPixels');
% pupilSize = Datapixx('GetPupilSizeSimple');
% [ppLeftMajor ppLeftMinor ppRightMajor ppRightMinor] = Datapixx('GetPupilSize');
% [ppLeftX ppLeftY ppRightX ppRightY] = Datapixx('GetPupilCoordinatesInPixels');
% [CRLeftX CRLeftY CRRightX CRRightY] = Datapixx('GetCRCoordinatesInPixels');
% Datapixx('SetupTPxSchedule', [bufferbaseAddress=12e6, numberOfEyeData=60000]);
% Datapixx('StopTPxSchedule');
% Datapixx('StartTPxSchedule');
% [bufferData, underflow, overflow] = Datapixx('ReadTPxData', numFrames);
% status = Datapixx('GetTPxStatus');
% Datapixx('ShowOverlay');
% Datapixx('HideOverlay');
% Datapixx('SetTPxSleep');
% Datapixx('SetTPxAwake');
% Datapixx('EnableSearchLimits');
% Datapixx('DisableSearchLimits');
% Datapixx('ClearSearchLimits');
% Datapixx('SetSearchLimits', leftEye, rightEye);
% [leftEye, rightEye] = Datapixx('GetSearchLimits');
% DEPRECATED: Datapixx('EnableTrackpixxAnalogOutput' ,[eyeNumber=0]);
% Datapixx('EnableTPxAnalogOut', [DAC0=0, DAC1=0, DAC2=0, DAC3=0]);
% DEPRECATED: Datapixx('DisableTrackpixxAnalogOutput');
% Datapixx('DisableTPxAnalogOut');
% Datapixx('PpSizeCalGetData');
% Datapixx('PpSizeCalGetDataComplete');
% Datapixx('PpSizeCalLinearRegression');
% Datapixx('PpSizeCalClear');
% Datapixx('PpSizeCalSet');
% [slope_r_x, slope_r_y, slope_l_x, slope_l_y] = Datapixx('PpSizeCalGet');
% fov_h = Datapixx('GetHorizontalFOV');
% fov_v = Datapixx('GetVerticalFOV');
% pxSize = Datapixx('GetPixelSize');
% pxDensity = Datapixx('GetPixelDensity');
% Datapixx('SetLens', lens);
% lens = Datapixx('GetLens');
% Datapixx('SetDistance', distance);
% distance = Datapixx('GetDistance');
% [leftFixationFlag, rightFixationFlag] = Datapixx('IsSubjectFixating');
% [leftSaccadeFlag, rightSaccadeFlag] = Datapixx('IsSubjectMakingSaccade');
% Datapixx('SetFixationThresholds' [, maxSpeed=2500] [, minNumberOfConsecutiveSamples=25]);
% Datapixx('SetSaccadeThresholds' [, minSpeed=10000] [, minNumberOfConsecutiveSamples=10]);
% 
% % TRACKPixx /mini Only Functions:
% time = Datapixx('GetTimeTPxMini');
% Datapixx('OpenTPxMini', [screenProportion=80]);
% Datapixx('CloseTPxMini');
% Datapixx('CalibrateTargetTPxMini', targetId);
% Datapixx('FinishCalibrationTPxMini');
% Datapixx('GetEyeImageTPxMini', imageArray);
% targets, targetsCount = Datapixx('InitializeCalibrationTPxMini');
% [bufferData] = Datapixx('ReadTPxMiniData');
% Datapixx('SetScreenProportionTPxMini', proportion]);
% Datapixx('SetupTPxMiniRecording', numberOfFrame);
% status = Datapixx('GetTPxMiniStatus');
% Datapixx('RecordTPxMiniData');
% Datapixx('StopTPxMiniRecording');
% [bufferData] = Datapixx('GetTPxMiniRecordedData');
% [bufferData] = Datapixx('GetTPxMiniLastEyePosition');
% 
% % Reading and writing register cache:
% Datapixx('RegWr');
% Datapixx('RegWrRd');
% Datapixx('RegWrVideoSync');
% Datapixx('RegWrRdVideoSync');
% Datapixx('RegWrPixelSync', pixelSequence [, timeout=255]);
% isTimeout = Datapixx('RegWrRdPixelSync, pixelSequence [, timeout=255]);
% 
% % Miscellaneous Routines:
% Datapixx('StopAllSchedules');
% error = Datapixx('GetError');
% Datapixx('ClearError');
% Datapixx('Reset');
% Datapixx('ResetAll');
% VERSION: 3.8.8370 [13/OCT/2020]
% Uses VPixx Device Server







