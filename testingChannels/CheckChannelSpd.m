% CheckChannelSpd
%
% This is to check channel power measuremtns.
%
% See Also:
%    MeasureChannelSpd.

% History:
%    12/01/21 smo       - Started on it for checking the spds of the new projector
%                         for SACC project.
%    03/28/22 smo       - Added an option to select the measurement device either
%                         PR670 or power meter.
%    03/29/22 dhb, smo  - Add in analysis.

%% Initialize.
clear; close all;

%% Set parameters.
S = [380 2 201];
nPrimaries = 3;
nChannels = 16;

targetChannels = [2 4 6 8 10 12];
nTargetChannels = length(targetChannels);

% Power meter settings.
powerMeterWl = 550;
wls = SToWls(S);
powerMeterWlIndex = find(wls == powerMeterWl);

projectorModeNormal = false;
VERBOSE = true;

%% Load spectrum data here.
DEVICE = 'PR670';

% Make a string for save file name.
% It used to be normal / steady-on
switch projectorModeNormal
    case true
        projectorMode = 'Normal';
    case false
        projectorMode = 'steady-on';
end

% Load the data here.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    testFilename = GetMostRecentFileName(testFiledir,sprintf('SpdData_%s_%s',DEVICE,projectorMode));
    prData = load(testFilename);
else
    error('Cannot find data file');
end

for pp = 1:nPrimaries
    targetChPeakWls(:,pp) = FindPeakSpds(prData.spdMeasured{pp},'verbose',false);
end

%% Load powermeter data here.
curDir = pwd;
cd(testFiledir);

% DATASET1
% Data with fixed wavelength sensitivity (550 nm).
powerSingleNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalSingle');
powerSingleSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnSingle');
powerWhiteNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalWhite');
powerWhiteSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnWhite');

% DATASET2
% Data with different wavelength.
% Make a string for save file name.
switch projectorModeNormal
    case true
        projectorMode = 'NormalMode';
    case false
        projectorMode = 'SteadyOnMode';
end

DEVICE = 'PowerMeter';
date = '0329';
fileType = '.csv';
targetChPeakWls = targetChPeakWls(:,3);
dataRange = 'D17:D17';

% Load sinlge peak data here.
for pp = 1:nPrimaries
    for cc = 1:nTargetChannels
        targetChPeakWl = targetChPeakWls(cc);
        targetCh = targetChannels(cc);
        fileName = append(DEVICE,'_',projectorMode,'_Primary',...
            num2str(pp),'_Ch',num2str(targetCh),'_',num2str(targetChPeakWl),'nm_',date,fileType);
        readFile = readmatrix(fileName, 'Range', dataRange);
        powerMeterWatt(cc,pp) = readFile;
    end
end

% Load white data here.
for cc = 1:nTargetChannels
    targetChPeakWl = targetChPeakWls(cc);
    fileName = append(DEVICE,'_',projectorMode,'_White_',num2str(targetChPeakWl),'nm_',date,fileType);
    readFile = readmatrix(fileName, 'Range', dataRange);
    powerMeterWhiteWatt(cc,:) = readFile;
end

%% Load power meter spectral sensitivity here.
powerMeterSensitivity = xlsread('PowerMeterResponsivityLocal.xlsx');
T_powerMeterRaw = SplineCmf(powerMeterSensitivity(:,1),powerMeterSensitivity(:,2)',S);

% Normalize power meter sensitivity so that wl set during measurement
% is one.
T_powerMeterMatch = T_powerMeterRaw/T_powerMeterRaw(powerMeterWlIndex);

% Plot it.
figure; clf; hold on;
plot(powerMeterWl,T_powerMeterMatch(powerMeterWlIndex),...
    'o','MarkerFaceColor','r','MarkerEdgeColor',zeros(3,1),'MarkerSize',7);
plot(wls,T_powerMeterMatch','k');
xlabel('Wavelength (nm)','FontSize',15');
ylabel('Normalized sensitivity', 'FontSize', 15);
legend('Normalized criteria');

% Back to the current path.
cd(curDir);

%% Power comes in as watts! Convert to mW.
wattToMWatt = 1000;

% DATASET1
powerSingleNormalMW = powerSingleNormalWatt .* wattToMWatt;
powerSingleSteadyOnMW = powerSingleSteadyOnWatt .* wattToMWatt;
powerWhiteNormalMW = powerWhiteNormalWatt .* wattToMWatt;
powerWhiteSteadyOnMW = powerWhiteSteadyOnWatt .* wattToMWatt;

% DATASET2
powerMeterMW = powerMeterWatt .* wattToMWatt;
powerMeterWhiteMW = powerMeterWhiteWatt .* wattToMWatt;

% Resize power meter data.
if (projectorModeNormal)
    powerMeterSingle = powerSingleNormalMW;
    powerMeterWhite = powerWhiteNormalMW;
else
    powerMeterSingle = powerSingleSteadyOnMW;
    powerMeterWhite = powerWhiteSteadyOnMW;
end

% Single peak data. Matching the shpae with the spectrum data.
powerMeterSingle = reshape(powerMeterSingle,length(targetChannels),nPrimaries);

%% Plot measured spectra.
if (VERBOSE)
    % Single peak spectrum.
    for pp = 1:nPrimaries
        figure; clf;
        plot(SToWls(S),prData.spdMeasured{pp});
        title(append('Screen Primary: ',num2str(pp),' ',projectorMode),'FontSize',15);
        xlabel('Wavelength (nm)','FontSize',15);
        ylabel('Spectral Irradiance','FontSize',15);
    end
    
    % White.
    figure; clf;
    plot(SToWls(S),prData.spdMeasuredWhite);
    title(append('White ',projectorMode),'FontSize',15);
    xlabel('Wavelength (nm)','FontSize',15);
    ylabel('Spectral Irradiance','FontSize',15);
end

%% Find scale factors for each measurement
for pp = 1:nPrimaries
    % Grab single primary spds
    spdSingle = prData.spdMeasured{pp};
    
    % Integrate against power meter sensitivity.
    F(:,pp) = (T_powerMeterMatch * spdSingle)';
end

FWhite = (T_powerMeterMatch * prData.spdMeasuredWhite);

% Compute factors k.
% These should be the same for all measurements.
k = powerMeterMW./F;
kWhite = powerMeterWhite/FWhite;

% Plot it.
figure; clf; hold on;
plot(kWhite,'bo','MarkerSize',12,'MarkerFaceColor','b');
plot(k,'ro','MarkerSize',12,'MarkerFaceColor','r');
xlabel('Target Channels','FontSize',15);
ylabel('Coefficient k','FontSize',15);
ylim([0 1.2*max(k(:))]);
legend('White','Single peak','FontSize',13);
