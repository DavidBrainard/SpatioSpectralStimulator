% QuantizationTest
%
% Work through exactly how things get quantized, so as to know enough to
% ensure consistency between physical display and our calculations.

% History:
%    11/05/21  dhb  Wrote it.

%% Clear
clear; close all;

%% Load in calibration file to play with
%
% This is a standard calibration file for the DLP projector,
% with the subprimaries set to something.  As we'll see below,
% we're going to rewrite those.nPrimaries
projectorCalName = 'SACC';
projectorNInputLevels = 256;

%% Load projector calibration.
projectorCal = LoadCalFile(projectorCalName);
screenCalObj = ObjectToHandleCalOrCalStruct(projectorCal);
CalibrateFitGamma(screenCalObj, projectorNInputLevels);
P_device = screenCalObj.get('P_device');
nPrimaries = size(P_device,2);

%% Refit gamma
gammaMethod = 'identity';
screenCalObj.set('gamma.fitType',gammaMethod);
CalibrateFitGamma(screenCalObj, projectorNInputLevels);

%% Set projector gamma method
%
% If we set to 0, there is no quantization and the result is excellent.
% If we set to 2, this is quantized at 256 levels and the result is more
% of a mess.  The choice of 2 represents what we think will actually happen
% since the real device is quantized.
%
% The point cloud method below reduces this problem.
projectorGammaMethod = 2;
SetGammaMethod(screenCalObj,projectorGammaMethod);

%% Set up finely spaced primary values and convert to settings.
%
% Find out what settings values come back
finePrimariesValues = linspace(0,1,2^16);
finePrimaries = finePrimariesValues(ones(nPrimaries,1),:);
fineSettings = PrimaryToSettings(screenCalObj,finePrimaries);
for ii = 1:size(fineSettings,2)
    if (any(fineSettings(:,ii) ~= fineSettings(1,ii)))
        error('Settings ');
    end
end

% Quick plot of what comes back for first primary channel.
% Look at expanded view in the plot
figure; clf;
subplot(1,2,1); hold on
plot(finePrimariesValues,fineSettings(1,:),'r+','MarkerSize',8);
xlabel('Primary Value');
ylabel('Setting Value');
axis('square'); xlim([0 0.1]); ylim([0 0.1]);
subplot(1,2,2); hold on
plot(finePrimariesValues,fineSettings(1,:),'r+','MarkerSize',8);
xlabel('Primary Value');
ylabel('Setting Value');
axis('square'); xlim([0.9 1]); ylim([0.9 1]);

% Check number of unique settings and what they are
expectedSettings = linspace(0,1,projectorNInputLevels);
for pp = 1:nPrimaries
    % Get and check number of unique settings
    uniqueSettingsValues(pp,:) = unique(fineSettings(pp,:));
    if (length(uniqueSettingsValues(pp,:)) ~= projectorNInputLevels)
        error('Do not get right number of unique setings values');
    end   
    if (any(uniqueSettingsValues(pp,:) ~= expectedSettings))
        error('Unique settings are not as expected');
    end
end

% Convert to integers
uniqueIntegerValues = SettingsToIntegers(uniqueSettingsValues,'nInputLevels',projectorNInputLevels);
temp1 = 0:(projectorNInputLevels-1);
temp2 = uniqueIntegerValues-temp1(ones(nPrimaries,1),:);
if (any(temp2(:) ~= 0))
    error('Conversion to integers not working as expected');
end

% If this conversion is all as expected, we shouldn't really have to round
% when we compute integers.  Let's check.
uniqueIntegerValuesNoRound = uniqueSettingsValues*(projectorNInputLevels-1);
roundCheck = max(abs(uniqueIntegerValuesNoRound(:)-round(uniqueIntegerValuesNoRound(:))));
if (roundCheck > 1e-10)
    error('Settings are not quantized the way we think they should be (first check)');
end
if (max(abs(uniqueIntegerValuesNoRound(:)-uniqueIntegerValues(:))) > 1e-10)
    error('Settings are not quantized the way we think they should be (second check)');
end
    
% Convert back to settings
settingsCheck = IntegersToSettings(uniqueIntegerValues,'nInputLevels',projectorNInputLevels);
if (max(abs(settingsCheck(:)-uniqueSettingsValues(:))) > 0)
    error('Cannot get settings back from integers');
end

