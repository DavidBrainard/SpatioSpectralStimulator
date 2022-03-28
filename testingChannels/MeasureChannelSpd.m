% MeasureChannelSpd
%
% This is to measure channel spd using either spectroradiometer or power
% meter.

% History:
%    12/01/21 smo  Started on it for checking the spds of the new projector
%                  for SACC project.
%    03/28/22 smo  Added an option to select the measurement device either
%                  PR670 or power meter.

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
VERBOSE = false;

%% Open the projector screen and measurement device ready.
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


%% Plot it.
if (VERBOSE)
    figure; clf; hold on;
    for pp = nPrimaries
        screenPrimaryStruct = append('screenPrimary',num2str(pp));
        subplot(nPrimaries,1,pp);
        
        for cc = 1:nMeasureChannels
            targetChannel = targetChannels(cc);
            plot(SToWls(S),spdMeasured(:));
            title('Raw');
            xlabel('Wavelength (nm)');
            ylabel('Spectral power distribution');
        end
    end
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
    % Make a string for save file name.
    switch projectorModeNormal
        case true
            projectorMode = 'Normal';
        case false
            projectorMode = 'SteadyOn';
    end
    
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

%% Save powermeter data here (temporary).
powerNormalRaw = xlsread('PowerMeterProcessedData.xlsx','NormalSingle');
powerSteadyRaw = xlsread('PowerMeterProcessedData.xlsx','SteadyOnSingle');

powerNormal = powerNormalRaw(:,3);
powerNormalSort = powerNormal(2:4:end);

powerSteady = powerSteadyRaw(:,3);

nDataPoints = 18;
powerNormalSort = 
