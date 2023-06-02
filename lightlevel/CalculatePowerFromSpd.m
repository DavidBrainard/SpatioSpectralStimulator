% CalculatePowerFromSpd.
%
% This calculates the luminance (cd/m2) from the spectrum.

% History:
%    6/2/23   smo     - Wrote it.

%% Initialize.
clear; close all;

%% Load the target spectrum.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'CheckCalibration');
    testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck');
    data = load(testFilename);
    spds = data.ptCldScreenSpdMeasuredCheckCal;
else
    error('Cannot find data file');
end

% Get spectrum range.
S = data.theData.S;

% Here we set the spectrum as the background of the stimuli (null
% stimulus).
spd = spds(:,1);

%% Load the power meter sensitivity.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    PWsensitivityFilename = fullfile(testFiledir,'PowerMeterResponsivityLocal.xlsx');
    powerMeterSensitivity = xlsread(PWsensitivityFilename);
else
    error('Cannot find data file');
end

% Match the wavelength range.
T_powerMeterRaw = SplineCmf(powerMeterSensitivity(:,1),powerMeterSensitivity(:,2)',S);

% Normalize power meter sensitivity accoridng to target wavelength when
% measuring the power meter. Default to 550 nm.
targetWl = 550;
wls = SToWls(S);
powerMeterWlIndex = find(wls == targetWl);
T_powerMeterMatch = T_powerMeterRaw/T_powerMeterRaw(powerMeterWlIndex);

%% Set coefficient value k.
%
% This is based on the pupil size of 4 mm which matches the situation of
% our experiment.
k = 0.0100;

%% Calculate luminance here.
%
% Power is in mW unit (milli watts).
power_mW = k * T_powerMeterMatch * spd;
fprintf('Calculated power is (%.6f) \n',power_mW);
