function [k] = SpdToPower(spd,powerWatt,options)
% This calculates the coefficient k that converts from relative irradiance
% (spectroradiometer) to the absolute (power meter).
%
% Syntax:
%    [k] = SpdToPower(spds,powerWatt)
%
% Description:
%    This calculates k which is the coefficient converting between the
%    relative irradiance and the absolute value.
%
% Inputs:
%    spd                        - dd
%    powerWatt                  -
%
%
% Optional key/value pairs:
%    S                          -
%    targetWls                  -
%    wattToMWatt                -
%    verbose                    - Default to true. Controls printout.
%
% See also:
%    MeasureChannelSpd, CheckChannelSpd.

% History:
%    03/31/22 smo               - Started on it.

%% Set parameters.
arguments
    spd
    powerWatt
    options.PWsensitivityFilename = 'PowerMeterResponsivityLocal.xlsx'
    options.S (1,3) = [380 2 201]
    options.targetWls {mustBeInRange(options.targetWls,380,780,"inclusive")} = 550
    options.wattToMWatt = 1000
    options.verbose (1,1) = false
end

%% Check the sizes of the input array.
nSpds = size(spd,2);
nPowerMeterMeasures = length(powerWatt);
if (~(nSpds == nPowerMeterMeasures))
    error('The size does not match between spds and power meter measurements!');
end

%% Load the power meter sensitivity.
%
% Load the raw sensitivity curve data.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    powerMeterSensitivity = xlsread(options.PWsensitivityFilename);
else
    error('Cannot find data file');
end

% Match the wavelength range.
T_powerMeterRaw = SplineCmf(powerMeterSensitivity(:,1),powerMeterSensitivity(:,2)',options.S);

% Normalize power meter sensitivity accoridng to target wavelength when
% measuring the power meter. Default to 550 nm.
nTargetWls = length(options.targetWls);
wls = SToWls(options.S);

powerMeterWlIndex = find(wls == options.targetWls);
T_powerMeterMatch = T_powerMeterRaw/T_powerMeterRaw(powerMeterWlIndex);

% Plot it.
if (options.verbose)
    figure; clf; hold on;
    plot(options.targetWls,T_powerMeterMatch(powerMeterWlIndex),...
        'o','MarkerFaceColor','r','MarkerEdgeColor',zeros(3,1),'MarkerSize',7);
    plot(wls,T_powerMeterMatch','k');
    xlabel('Wavelength (nm)','FontSize',15');
    ylabel('Normalized sensitivity', 'FontSize', 15);
    legend('Normalized criteria');
end

%% Integrate against power meter sensitivity.
%
% Scale power meter measurements in milli watt unit.
powerMWatt = powerWatt .* options.wattToMWatt;

% Find scale factors for each measurement.
%
% Here we calculate F which is basically multiplying power meter
% sensitivity and spds.
F = T_powerMeterMatch * spd;

%% Compute factors k.
%
% The factor k is the coefficient that converts F value calculated above to
% power meter measurements in Watt.
k = powerMWatt./F;

end
