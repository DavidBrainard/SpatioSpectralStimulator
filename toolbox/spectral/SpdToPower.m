function [k] = SpdToPower(spd,powerWatt,options)
% This calculates the coefficient k that converts from relative irradiance
% (spectroradiometer) to the absolute (power meter).
%
% Syntax:
%    [k] = SpdToPower(spd,powerWatt)
%
% Description:
%    This calculates k which is the coefficient converting between the
%    relative irradiance and the absolute value. 
%
%    It calculates one k value at a time, so if you want to calculate for
%    multiple cases, you need to make a loop in your run script.
%
% Inputs:
%    spd                        - Target spd measured by spectroradiometer.
%    powerWatt                  - Target power measured using power meter
%                                 in Watt units. 
%
% Outputs:
%    k                          - The conversion factor.  This converts
%                                 the raw radiometer values to power in
%                                 units defined by the wattToMWatt
%                                 key/value pair.  By default, this has
%                                 value 1000, so that k converts measured
%                                 spectra to power in mWatt units.  Note
%                                 that spectra here follow the PTB
%                                 convention of power per wavelength band.
%
% Optional key/value pairs:
%    S                          - Default to [380 2 201]. Wavelength range.
%    targetWl                   - Default to 550. The sensitivity of power
%                                 meter set for the measurement. 
%    wattToMWatt                - Default to 1000. It converts the power
%                                 unit from Watt to MilliWatt. 
%    verbose                    - Default to false. Controls printout.
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
    options.targetWl {mustBeInRange(options.targetWl,380,780,"inclusive")} = 550
    options.wattToMWatt = 1000
    options.verbose (1,1) = false
end

%% Check the sizes of the input array.
nSpd = size(spd,2);
nPowerMeasure = length(powerWatt);
if (~(nSpd == nPowerMeasure))
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
wls = SToWls(options.S);
powerMeterWlIndex = find(wls == options.targetWl);
T_powerMeterMatch = T_powerMeterRaw/T_powerMeterRaw(powerMeterWlIndex);

% Plot it.
if (options.verbose)
    figure; clf; hold on;
    plot(options.targetWl,T_powerMeterMatch(powerMeterWlIndex),...
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
