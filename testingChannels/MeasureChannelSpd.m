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

%% Initialize.
clear; close all;

%% Set parameters.
S = [380 2 201];
nPrimaries = 3;
nChannels = 16;
screenSettings = [1 1 1]';

% Channel settings.
targetChannels = [2 4 6 8 10 12];
targetChannelSettingValue = 1;
logicalToPhysical = [0:15];

% Black correction for spectrum measurements.
arbitraryBlack = 0.05;

% Choose which device to use within [PR670, powermeter].
DEVICE = 'PR670';
projectorModeNormal = false;
powerMeterWaitTimeSec = 5;
VERBOSE = true;

% Make a string for save file name.
switch projectorModeNormal
    case true
        projectorMode = 'Normal';
    case false
        projectorMode = 'SteadyOn';
end

%% Open the projector screen and measurement device ready.
[window, windowRect] = OpenPlainScreen(screenSettings,...
    'projectorMode',projectorModeNormal);

if (strcmp(DEVICE,'PR670'))
    OpenSpectroradiometer;
end

%% Get ready.
initialDealySec = 5;
fprintf('Measurement will begin in (%d) seconds...\n',initialDealySec);
WaitSecs(initialDealySec);

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

%% This part measures arbitrary black.
switch DEVICE
    case 'PR670'
        if (~arbitraryBlack == 0)
            channelSettingsBlack = ones(nChannels,nPrimaries) .* arbitraryBlack;
            SetChannelSettings(channelSettingsBlack);
            GetChannelSettings;
            spdMeasuredBlack = MeasureSpectroradiometer;
            disp('Black measurement has been complete!');
        else
            % If the black level is not set, just pass the zero spectrum.
            spdMeasuredBlack = zeros(S(3),1);
        end
        
    case 'powermeter'
        % Set black back to zero. We don't use arbitrary black for power
        % meter measurements. Instead, power meter is calibrated to the
        % screen with all LED lights off.
        arbitraryBlack = 0.0;
    otherwise
        error('Device should be selected either PR670 or powermeter');
end

%% This part measure singe channels.
%
% Set channel settings and measure here. All data will be stored in
% spdmeasured.
nMeasureChannels = length(targetChannels);

for pp = 1:nPrimaries
    for cc = 1:nMeasureChannels
        % Set channel settings here.
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
                WaitSecs(powerMeterWaitTimeSec);
            otherwise
        end
    end
    
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
