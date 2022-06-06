% SimulateScreenPrimaries
%
% This is to simulate a spectrum using the LED channels. This is useful to
% make some spectrum graphs for graphical presentations.
%
% Current version can generate spectrum that has the same CIE XYZ values.
% The result is depending on how to set the primaries of the test display.

% History
%    6/5/22    smo    - Written it.

%% Initialize.
clear; close all;

%% Load the spectum data here.
calFilename = 'SACCPrimary3';
if (ispref('BrainardLabToolbox','CalDataFolder'))
    testFiledir = getpref('BrainardLabToolbox','CalDataFolder');
    allCalData = load(fullfile(testFiledir,calFilename));
else
    error('Cannot find data file');
end

% Get the most recent data file.
mostRecentIndex = size(allCalData.cals,2);
calData = allCalData.cals{mostRecentIndex};

% Read out some variables.
B_channels = calData.processedData.P_device;
S = calData.rawData.S;
wls = SToWls(S);

% Plot all channels.
figure; hold on;
plot(wls,B_channels);
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
xlim([380 780]);

%% Load CMFs here.
%
% Load cone sensitivity function.
load T_cones_sp
T_Cones = T_cones_sp;
T_Cones = SplineCmf(SToWls(S_cones_sp),T_Cones,wls)';

% Load color matching function.
load T_xyzJuddVos
T_XYZ = T_xyzJuddVos;
T_XYZ = SplineCmf(SToWls(S_xyzJuddVos),T_XYZ,wls)';

% Conventional CRT Monitor.
load B_monitor
B_monitor = SplineSpd(SToWls(S_monitor),B_monitor,wls);

%% Make a test display with modulated primaries.
nChannels = size(B_channels,2);
testDisplayNum = 1;
switch testDisplayNum
    case 1
        primary1Channels = [12 13 15 16]; % No 14
        primary2Channels = [6 7 8 9 10];
        primary3Channels = [1 2 3];
    case 2
        primary1Channels = [12]; % No 14
        primary2Channels = [9];
        primary3Channels = [2];
    case 3
        primary1Channels = [8 9 10 11 12 13 15 16]; % No 14
        primary2Channels = [4 5 6 7];
        primary3Channels = [1 2 3];
    otherwise
end

% Set test display primaries.
B_primary1 = sum(B_channels(:,primary1Channels),2);
B_primary2 = sum(B_channels(:,primary2Channels),2);
B_primary3 = sum(B_channels(:,primary3Channels),2);
B_testDisplay = [B_primary1 B_primary2 B_primary3];

% Plot the modulated primaries.
figure; hold on;
plot(wls,B_testDisplay(:,1),'r-');
plot(wls,B_testDisplay(:,2),'g-');
plot(wls,B_testDisplay(:,3),'b-');
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
xlim([380 780]);
legend('Primary1','Primary2','Primary3');
% title('Test display primaries');

%% Find a desired spectrum.
%
% Set target XYZ values here.
targetXYZ = [100 100 100];

% Find the linear transformation M_XYZtoCones that maps XYZ tristimulus
% coordinates to cone responses.
M_XYZtoCones = (T_XYZ\T_Cones);

% Find dRGB input ratio for monitor to simulate D65
M_XYZtoRGB = targetXYZ * inv(B_testDisplay' * T_XYZ);
targetSpd = M_XYZtoRGB * B_testDisplay';

% Check if XYZ values were calculated back okay.
checkXYZ = targetSpd * T_XYZ;

% Plot the spectrum of desired XYZ values.
figure;
plot(wls,targetSpd);
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
xlim([380 780]);
ylim([0 8]);
title('A spectrum found matching XYZ');

%% ipRGCs example.
numChannelPeak = 3;
offChannels = [4,7:12,14:16]
onChannels  = [offChannels numChannelPeak]

B_offPeak = sum(B_channels(:,offChannels),2)
B_onPeak = sum(B_channels(:,onChannels),2)

figure; hold on;
plot(wls,B_offPeak,'k-','LineWidth',1)
plot(wls,B_onPeak,'g--','LineWidth',1)
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
legend('Off-peak','On-peak')

%% Addtivity example.
%
% Set which channels to test.
additivityChannels = [2 7 13];

% Make spectra here.
B_addtivity = B_channels(:,additivityChannels);
B_addtivitySum = sum(B_addtivity,2);

% Plot separate single peaks.
figure;
plot(wls,B_addtivity)
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
xlim([380 780]);

% Plot sum of single peaks.
figure;
plot(wls,B_addtivitySum);
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');
xlim([380 780]);
