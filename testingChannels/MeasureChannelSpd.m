% MeasureChannelSpd
%
% This is to measure channel spd using either spectroradiometer or power
% meter.

% History:
%    12/01/21 smo  Started on it for checking the spds of the new projector
%                  for SACC project.
%    03/28/22 smo  Added an option to select the measurement device either
%                  PR670 or power meter.

%% Initialize.
clear;  close all;

%% Set parameters.
S = [380 2 201];
nPrimaries = 3;
nChannels = 16;

targetChannels = [2 4 6 8 10 12];
targetChannelSettingValue = 1;
logicalToPhysical = [0:15];

screenSettings = [1 1 1]';

% Choose which device to use within [PR670, powerMeter].
DEVICE = 'PR670';
projectorModeNormal = false;
powerMeterWaitTimeSec = 5;
VERBOSE = true;

% Make a string for save file name.
switch projectorModeNormal
    case true
        projectorMode = 'normal';
    case false
        projectorMode = 'steady-on';
end

LOADDATA = true;

%% Load the data if you want.
if (LOADDATA)
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('SpdData_%s_%s',DEVICE,projectorMode));
        load(testFilename);
    else
        error('Cannot find data file');
    end
end

%% Open the projector screen and measurement device ready.
if (~LOADDATA)
    [window, windowRect] = OpenPlainScreen(screenSettings,...
        'projectorMode',projectorModeNormal);
    
    if (strcmp(DEVICE,'PR670'))
        OpenSpectroradiometer;
    end
    
    %% This part measures full white.
    channelSettingsWhite = ones(nChannels,nPrimaries);
    SetChannelSettings(channelSettingsWhite);
    GetChannelSettings;
    
    % Measure it.
    switch DEVICE
        case 'PR670'
            spdMeasuredWhite = MeasureSpectroradiometer;
        case 'powermeter'
            WaitSecs(powerMeterWaitTimeSec);
        otherwise
    end
    disp('White measurement has been complete!');
    
    %% This part measure singe channels.
    %
    % Set channel settings and measure here. All data will be stored in
    % spdmeasured.
    nMeasureChannels = length(targetChannels);
    
    for pp = 1:nPrimaries
        for cc = 1:nMeasureChannels
            % Set channel settings here.
            targetChannel = targetChannels(cc);
            channelSettings = zeros(nChannels,nPrimaries);
            channelSettings(targetChannel, pp) = targetChannelSettingValue;
            SetChannelSettings(channelSettings);
            
            % Verbose.
            fprintf('Screen primary (%d) / channel (%d) is now displaying \n', pp, targetChannel);
            
            % Check the current channel settings.
            GetChannelSettings;
            
            % Measure it.
            switch DEVICE
                case 'PR670'
                    spdChannelSingles(:,cc) = MeasureSpectroradiometer;
                case 'powermeter'
                    WaitSecs(powerMeterWaitTimeSec);
                otherwise
            end
        end
        spdMeasured{pp} = spdChannelSingles;
    end
    
    %% Close.
    if (strcmp(DEVICE,'powermeter'))
        CloseScreen;
    end
    
    if (strcmp(DEVICE,'PR670'))
        CloseSpectroradiometer;
    end
    
    %% Save the data.
    if (strcmp(DEVICE,'PR670'))
        % Save the file here.
        if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
            testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
            dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            testFilename = fullfile(testFiledir,sprintf('SpdData_%s_%s_%s',DEVICE,projectorMode,dayTimestr));
            save(testFilename,'S','spdMeasured','spdMeasuredWhite','targetChannels','targetChannelSettingValue');
            disp('Data has been saved successfully!');
        end
    end
    
    %% Say goodbye.
    disp('All measurement has been complete!');
end

%% Plot it.
if (VERBOSE)
    % Single peak spectrum.
    for pp = 1:nPrimaries
        figure; clf;
        plot(SToWls(S),spdMeasured{pp});
        title(append('Screen Primary: ',num2str(pp),' ',projectorMode),'FontSize',15);
        xlabel('Wavelength (nm)','FontSize',15);
        ylabel('Spectral Irradiance','FontSize',15);
    end
    
    % White.
    figure; clf;
    plot(SToWls(S),spdMeasuredWhite);
    title(append('White ',projectorMode),'FontSize',15);
    xlabel('Wavelength (nm)','FontSize',15);
    ylabel('Spectral Irradiance','FontSize',15);
end

%% Load powermeter data here (temporary).
powerSingleNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalSingle');
powerSingleSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnSingle');
powerWhiteNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalWhite');
powerWhiteSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnWhite');

wattToMWatt = 1000;

powerSingleNormalMW = powerSingleNormalWatt .* wattToMWatt;
powerSingleSteadyOnMW = powerSingleSteadyOnWatt .* wattToMWatt;
powerWhiteNormalMW = powerWhiteNormalWatt .* wattToMWatt;
powerWhiteSteadyOnMW = powerWhiteSteadyOnWatt .* wattToMWatt;

% Plot it.
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

%% Take the irradiance at target wavelength.
targetWls = 550;
Wls = SToWls(S);
targetWlsIndex = find(Wls == targetWls);

for pp = 1:nPrimaries
    spdSingle = spdMeasured{pp};
    spdSingleTargetWls(:,pp) = spdSingle(targetWlsIndex,:);
end

spdSingleTargetWls = reshape(spdSingleTargetWls, size(spdSingleTargetWls,1)*size(spdSingleTargetWls,2), 1);