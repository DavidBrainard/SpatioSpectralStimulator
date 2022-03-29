% CheckChannelSpd
%
% This is to check channel power measuremtns.

% History:
%    12/01/21 smo  Started on it for checking the spds of the new projector
%                  for SACC project.
%    03/28/22 smo  Added an option to select the measurement device either
%                  PR670 or power meter.
%    03/29/22 dhb, smo  Add in analysis.

%% Initialize.
clear;  close all;

%% Set parameters.
S = [380 2 201];
nPrimaries = 3;
nChannels = 16;

% Power meter
powerMeterWl = 550;
wls = SToWls(S);
powerMeterWlIndex = find(wls == powerMeterWl);

targetChannels = [2 4 6 8 10 12];
targetChannelSettingValue = 1;
logicalToPhysical = [0:15];

screenSettings = [1 1 1]';

% Choose which device to use within [PR670, powerMeter].
DEVICE = 'PR670';
projectorModeNormal = true;
powerMeterWaitTimeSec = 5;
VERBOSE = true;

% Make a string for save file name.
switch projectorModeNormal
    case true
        projectorMode = 'normal';
    case false
        projectorMode = 'steady-on';
end


%% Load the data if you want.
LOADDATA = true;
if (LOADDATA)
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('SpdData_%s_%s',DEVICE,projectorMode));
        prData = load(testFilename);
    else
        error('Cannot find data file');
    end

    %% Load powermeter data here (temporary).
    curDir = pwd;
    cd(testFiledir);
    powerSingleNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalSingle');
    powerSingleSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnSingle');
    powerWhiteNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalWhite');
    powerWhiteSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnWhite');

    % Read power meter spectral sensitivity
    powerMeterSensitivity = xlsread('PowerMeterResponsivityLocal.xlsx');
    T_powerMeter = SplineCmf(powerMeterSensitivity(:,1),powerMeterSensitivity(:,2)',S);

    % Normalize power meter sensitivity so that wl set during measurement
    % is one.
    T_powerMeter = T_powerMeter/T_powerMeter(powerMeterWlIndex);
    figure; clf; plot(wls,T_powerMeter','k');

    cd(curDir);
end

%% Power comes in as watts! Convert to mW.
wattToMWatt = 1000;
powerSingleNormalMW = powerSingleNormalWatt .* wattToMWatt;
powerSingleSteadyOnMW = powerSingleSteadyOnWatt .* wattToMWatt;
powerWhiteNormalMW = powerWhiteNormalWatt .* wattToMWatt;
powerWhiteSteadyOnMW = powerWhiteSteadyOnWatt .* wattToMWatt;

%% Plot measured spectra
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

%% Plot power meter data
if (VERBOSE)
    % Get the projector mode info.
    if (projectorModeNormal)
        powerMeterSingle = powerSingleNormalMW;
        powerMeterWhite = powerWhiteNormalMW;
    else
        powerMeterSingle = powerSingleSteadyOnMW;
        powerMeterWhite = powerWhiteSteadyOnMW;
    end
    
    % Single peak data.
    powerMeterSingle = reshape(powerMeterSingle,length(targetChannels),nPrimaries);
    
    for pp = 1:nPrimaries
        figure; clf;
        plot(powerMeterSingle(:,pp),'ro')
        xlabel('Target channel','FontSize',15);
        ylabel('Power meter (mw)','FontSize',15);
        title(append('Screen Primary: ',num2str(pp),' ',projectorMode),'FontSize',15);
    end
    
    % White data.
    figure; clf;
    plot(powerMeterWhite,'bo');
    xlabel('Target channel','FontSize',15);
    ylabel('Power meter (mw)','FontSize',15);
    title(append('White ',projectorMode),'FontSize',15);
end

% Find scale factors for each measurement
for pp = 1:nPrimaries
    % Grab single primary spds
    spdSingle = prData.spdMeasured{pp};

    % Integrate against power meter sensitivity
    F(:,pp) = (T_powerMeter*spdSingle)';
end

% Compute factors k.  These should be the same for all measurements
k = powerMeterSingle./F;

figure; clf; plot(k,'ro','MarkerSize',12,'MarkerFaceColor','r');
ylim([0 1.2*max(k(:))]);
