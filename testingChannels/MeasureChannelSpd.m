% MeasureChannelSpd
%
% This is to measure channel spd using either spectroradiometer or power
% meter.
%
% See Also:
%    CheckChannelSpd.

% History:
%    12/01/21 smo  - Started on it for checking the spds of the new projector
%                    for SACC project.
%    03/28/22 smo  - Added an option to select the measurement device either
%                    PR670 or power meter.
%    03/30/22 smo  - Deleted the option loading the data and added black
%                    correction for spectrum measurement.
%    04/19/22 smo  - Updates on setting arbitrary black.
%    05/05/22 smo  - Added an option to warm up the screen before starting
%                    the measurements.

%% Initialize.
clear; close all;

%% Set parameters.
S = [380 2 201];
nPrimaries = 3;
nChannels = 16;
screenSettings = [1 1 1]';
warmupTimeMin = 0;

% Channel settings.
targetChannels = [2 4 6 8 10 12];
targetChannelSettingValue = 1;
logicalToPhysical = [0:15];

% Black correction for spectrum measurements.
arbitraryBlack = 0.05;

% Choose which device to use within [PR670, powermeter].
DEVICE = 'powermeter';
projectorModeNormal = false;
powerMeterWaitTimeSec = 10;
VERBOSE = true;

% Make a string for save file name.
switch projectorModeNormal
    case true
        projectorMode = 'NormalMode';
    case false
        projectorMode = 'SteadyOnMode';
end

%% Warmup the screen if you want.
if (~warmupTimeMin == 0)
    WarmupScreen('projectorMode',projectorModeNormal,...
        'warmupTimeMin',warmupTimeMin,'verbose',VERBOSE);
end

%% Open the projector screen and measurement device ready.
%
% Open the projector with the set screen settings.
[window, windowRect] = OpenPlainScreen(screenSettings,...
    'projectorMode',projectorModeNormal);

% Get spectroradiometer ready.
if (strcmp(DEVICE,'PR670'))
    OpenSpectroradiometer;
end

%% Make a delay before starting measurements to get out of the room.
initialDelaySec = 3;
fprintf('Measurement will begin in (%d) seconds...\n',initialDelaySec);
if (strcmp(DEVICE,'powermeter'))
    disp('You should run the power meter program when the timer is done!');
end

% Display timer progress.
for tt = 1:initialDelaySec
    fprintf('Timer counting: (%d/%d) seconds...\n',tt,initialDelaySec);
    pause(1);
end

% Power meter turn on timer. You should turn on the power meter program as
% soon as this timer count starts.
if (strcmp(DEVICE,'powermeter'))
    disp('Now start power meter measurement on the program!');
    halfPowerMeterWaitTimeSec = powerMeterWaitTimeSec/2;
    for tt = 1:halfPowerMeterWaitTimeSec
        fprintf('Timer counting: (%d/%d) seconds...\n',tt,halfPowerMeterWaitTimeSec);
        pause(1);
    end
end

%% Measure full white.
channelSettingsWhite = ones(nChannels,nPrimaries);
SetChannelSettings(channelSettingsWhite);
GetChannelSettings;

% Measure it.
switch DEVICE
    case 'PR670'
        spdMeasuredWhite = MeasureSpectroradiometer;
    case 'powermeter'
        pause(powerMeterWaitTimeSec);
    otherwise
end
disp('White measurement has been complete!');

%% Measure single channel spectrum.
%
% All data will be stored in spdmeasured.
nMeasureChannels = length(targetChannels);

for pp = 1:nPrimaries
    for cc = 1:nMeasureChannels
        % Set channel settings here.
        if (strcmp(DEVICE,'powermeter'))
            clear arbitraryBlack;
            arbitraryBlack = 0;
        end
        targetChannel = targetChannels(cc);
        channelSettings = ones(nChannels,nPrimaries) .* arbitraryBlack;
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
                pause(powerMeterWaitTimeSec);
            otherwise
        end
        
        % This part measures arbitrary black for black correction.
        if (strcmp(DEVICE,'PR670'))
            if (~arbitraryBlack == 0)
                % Measure black level if the arbitrary black is set.
                channelTurnOff = 0;
                channelSettingsBlack = ones(nChannels,nPrimaries) .* arbitraryBlack;
                channelSettingsBlack(targetChannel,pp) = channelTurnOff;
                SetChannelSettings(channelSettingsBlack);
                GetChannelSettings;
                spdMeasuredBlack(:,cc) = MeasureSpectroradiometer;
                disp('Black measurement has been complete!');
            else
                % If the black level is not set, just pass the zero spectrum.
                spdMeasuredBlack = zeros(S(3),1);
            end
        end
    end
    
    % Black correction here.
    if(strcmp(DEVICE,'PR670'))
        spdMeasured{pp} = spdChannelSingles - spdMeasuredBlack;
    end
end

%% Find peak spectrum.
if (strcmp(DEVICE,'PR670'))
    for pp = 1:nPrimaries
        peakSpd(:,pp) = FindPeakSpds(spdMeasured{pp},'verbose',false);
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
    % Save the file here.
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,sprintf('SpdData_%s_%s_%s',DEVICE,projectorMode,dayTimestr));
        save(testFilename,'S','spdMeasured','spdMeasuredWhite','spdMeasuredBlack','targetChannels','targetChannelSettingValue','peakSpd');
        disp('Data has been saved successfully!');
    end
end

%% Say goodbye.
disp('All measurement has been complete!');

%% Plot it.
if (strcmp(DEVICE,'PR670'))
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
end
