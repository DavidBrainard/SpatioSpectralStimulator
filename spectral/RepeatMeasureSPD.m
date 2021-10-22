% RepeatMeasureSPD
% 
% This measures the SPD of the desired projector setting with repeatitions.
%
% 10/22/2021 smo Strated on it.
%
%% Initialize.
close all; clear all; clc;

%% Load spectral data you want to measure
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,'testImageData1');
    theData = load(testFilename);
end

%% Set parameters.
logicalToPhysical = [0:7 9:15];
subprimaryNInputLevels = 253;
nSubprimaries = 15;
nPrimaries = 3;
S = theData.S;
wls = SToWls(S);

%% Set projector and measurement device ready.
addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx'));
    
% Connect to the projector.
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

% Connect to PR-670
command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
CMCheckInit(5); 

% Set projector state.
command_primaries = 'vputil rw 0x1c8 0x0 -q quit'; % Set to 'Normal mode'
unix(command_primaries)

%% Set the desired spd to measure.
targetSpd = PrimaryToSpd(theData.subprimaryCalObjs{1},theData.projectorPrimaryPrimaries(:,1));
targetSubprimarySettings = zeros(nSubprimaries,nPrimaries);
targetSubprimarySettings(:,1) = PrimaryToSettings(theData.subprimaryCalObjs{1},theData.projectorPrimaryPrimaries(:,1));

% Set the projector subprimaries here.
for ss = 1:nSubprimaries
    Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(targetSubprimarySettings(ss,1)*(subprimaryNInputLevels-1))); % Primary 1
    Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(targetSubprimarySettings(ss,2)*(subprimaryNInputLevels-1))); % Primary 2
    Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(targetSubprimarySettings(ss,3)*(subprimaryNInputLevels-1))); % Primary 3
end

% Check the current projector currents levels
for j=1:3; % PrimaryColor (0-2)
    for i=1:16; % SubColor (0-15)
    currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
    end
end 
currents

%% Display plain screen on DLP using PTB.
PsychDefaultSetup(2); % PTB pre-setup
screens = Screen('Screens');
screenNumber = max(screens);
white = WhiteIndex(screenNumber);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

%% Measure it.
nMeas = 100; % Number of repetition for the measurements.
for i = 1:nMeas
   spdMeasured(:,i) = MeasSpd(S,5,'all');
   fprintf('           Measurement complete! - (%d/%d)\n',i,nMeas);
end

% Plot it.
figure; hold on;
plot(wls,spdMeasured);
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
