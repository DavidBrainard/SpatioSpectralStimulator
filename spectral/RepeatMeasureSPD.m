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

% Set projector state.
command_primaries = 'vputil rw 0x1c8 0x0 -q quit'; % Set to 'Normal mode'
unix(command_primaries)

% Connect to PR-670
command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
CMCheckInit(5); 

%% Set the desired spd to measure.
targetSpd = PrimaryToSpd(theData.subprimaryCalObjs{1},theData.projectorPrimaryPrimaries(:,1));
targetSubprimarySettings = zeros(nSubprimaries,nPrimaries);
for pp = 1:3
    targetSubprimarySettings(:,pp) = PrimaryToSettings(theData.subprimaryCalObjs{pp},theData.projectorPrimaryPrimaries(:,pp));
end

% % Check the current projector currents levels
% for j=1:3; % PrimaryColor (0-2)
%     for i=1:16; % SubColor (0-15)
%     currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
%     end
% end 
% currents

%% Display plain screen on DLP using PTB and measure it.
% Set number of repeatitions.
nMeas = 50; % Number of repetition for the measurements.

%% ver1
% for i = 1:nMeas
%     PsychDefaultSetup(2); % PTB pre-setup
%     screens = Screen('Screens');
%     screenNumber = max(screens);
%     white = WhiteIndex(screenNumber);
%     [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
% 
%     % Set the projector subprimaries here. All zero for start.
% for ss = 1:nSubprimaries
%     Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), 0); % Primary 1
%     Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), 0); % Primary 2
%     Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), 0); % Primary 3
% end   
%     % Set the projector subprimaries here.
% for ss = 1:nSubprimaries
%     Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(targetSubprimarySettings(ss,1)*(subprimaryNInputLevels-1))); % Primary 1
%     Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(targetSubprimarySettings(ss,2)*(subprimaryNInputLevels-1))); % Primary 2
%     Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(targetSubprimarySettings(ss,3)*(subprimaryNInputLevels-1))); % Primary 3
% end
%     spdMeasured(:,i) = MeasSpd(S,5,'all');
%     fprintf('           Measurement complete! - (%d/%d)\n',i,nMeas);
%     sca;
% end

%% ver2
for mm = 1:nMeas
    for pp = 1:nPrimaries
        % Display plain screen on DLP using PTB.
        PsychDefaultSetup(2); % PTB pre-setup
        screens = Screen('Screens');
        screenNumber = max(screens);
        white = WhiteIndex(screenNumber);
        [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
  
        otherPrimaries = setdiff(1:nPrimaries,pp);
    
    % Set projector current levels as the above settings.
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', pp-1, logicalToPhysical(ss), round(targetSubprimarySettings(ss,pp)*(subprimaryNInputLevels-1)));   % Target primary
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(1)-1, logicalToPhysical(ss), 0); % Other Primary 1
        Datapixx('SetPropixxHSLedCurrent', otherPrimaries(2)-1, logicalToPhysical(ss), 0); % Other Primary 2
    end
    
    % Check the current projector currents levels
    for j=1:3; % PrimaryColor (0-2)
        for i=1:16; % SubColor (0-15)
        currents(j,i) = Datapixx('GetPropixxHSLedCurrent',j-1,i-1);
        end
    end 
    currents

    % Measurement.
    spdMeasured(:,pp+nPrimaries*(mm-1)) = MeasSpd(S,5,'all');
    fprintf('           Measurement complete! - Primary %d\n',pp);
    
    % Close PTB screen.
    sca;
    
    end
    fprintf('           Measurement complete! - (%d/%d)\n',mm,nMeas);
end

%% Plot it.
figure; hold on;
plot(wls,spdMeasured);
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
