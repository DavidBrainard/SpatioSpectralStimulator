% MeasureChannelSpdQuick
%
% This is measure channel spd for the purpose of checking it quick.

% History:
%    12/01/21 smo  Started on it for checking the spds of the new projector
%                  for SACC project.

%% Open.
OpenPlainScreen([1 1 1]);
OpenSpectroradiometer;

% Set parameters.
S = [380 2 201];
channelSettingValue = 1;
nChannels = 16;
screenPrimary = 3;

% Measure it.
for cc = 1:nChannels;
    channelSettings = zeros(16,3);
    
    channelSettings(cc,screenPrimary) = channelSettingValue;
    SetChannelSettings(channelSettings);
    
    spdMeasured(:,cc) = MeasureSpectroradiometer;
end

for cc = 1:16
    spdNormalized(:,cc) = spdMeasured(:,cc) ./ max(spdMeasured(:,cc));
end

% Make a plot.
figure; clf;
subplot(2,1,1);
plot(SToWls(S),spdMeasured);
title('Raw');
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');

subplot(2,1,2);
plot(SToWls(S),spdNormalized);
title('Normalized');
xlabel('Wavelength (nm)');
ylabel('Spectral power distribution');

%% Close.
CloseScreen;
CloseSpectroradiometer;

%% Save the data.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    testFilename = fullfile(testFiledir,sprintf('channelSpdCheck_%s',dayTimestr));
    save(testFilename,'S','spdMeasured','spdNormalized','screenPrimary','channelSettingValue');
    disp('Data has been saved successfully!');
end
