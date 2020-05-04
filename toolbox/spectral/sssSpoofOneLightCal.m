function cal = sssSpoofOneLightCal(varargin)
% Make a OneLight calibration structure for the spatial-spectral simulator
%
% Syntax:
%    cal = sssSpoofOneLightCal
%
% Description:
%    Make a OneLight calibration structure that describes the spatial-spectral simulator.
%
% Inputs:
%    None.
%
% Outputs:
%    cal        - The generated calibration structure.
%
% Optional key/value pairs:
%    plotBasis  - Boolean (default false).  Make a plot of the basis.
%
% See Also:
%    sssGetSpectralBasis
%

% History:
%    05/03/20  dhb  Wrote it.


% Examples:
%{
    cal = sssSpoofOneLightCal;
%}

% Massage varargin and parse
p = inputParser;
%p.addParameter('plotBasis', false, @islogical);
p.parse(varargin{:});

%% ID
cal.calID = 'sssSpoofCal';

%% Describe
cal.describe.gammaFitType = 'betacdfpiecelin';
cal.describe.useAverageGamma = 1;
cal.describe.nGammaFitLevels = 256;
cal.describe.nGammaLevels = 24;
cal.describe.gammaNumberWlUseIndices = 5;
cal.describe.whichAverageGamma = 'median';
cal.describe.S = [380 2 201];
cal.describe.useOmni = false;

%% Spectral basis
B = sssGetSpectralBasis(cal.describe.S);
wls = SToWls(cal.describe.S);
cal.describe.numWavelengthBands = size(B,2);
cal.describe.nGammaBands = cal.describe.numWavelengthBands;
cal.raw.darkMeas = zeros(size(wls,1),2);
cal.raw.darkMeasCheck = cal.raw.darkMeas;
cal.raw.lightMeas = B;

%% Gamma
cal.describe.gamma.gammaLevels = linspace(1/cal.describe.nGammaLevels,1,cal.describe.nGammaLevels);
cal.describe.gamma.gammaBands = 1:cal.describe.nGammaBands;
for ii = 1:cal.describe.numWavelengthBands
    for jj = 1:cal.describe.nGammaLevels
        cal.raw.gamma.rad(ii).meas(:,jj) = cal.describe.gamma.gammaLevels(jj)*cal.raw.lightMeas(:,ii);
    end
end

%% Get cal.computed

% Wavelength sampling
cal.computed.pr650S = cal.describe.S;
cal.computed.pr650Wls = SToWls(cal.computed.pr650S);
cal.computed.commonWls = cal.computed.pr650Wls;

% Copy over spectral data
cal.computed.pr650M = cal.raw.lightMeas;
cal.computed.pr650Md = cal.raw.darkMeas;
cal.computed.pr650MeanDark = mean(cal.computed.pr650Md,2);
cal.computed.pr650MeanDark(cal.computed.pr650MeanDark < 0) = 0;

% Use the stops information to get the fraction of max power for each gamma measurement, in range [0,1].
cal.computed.gammaInputRaw = [0 ; cal.describe.gamma.gammaLevels'];

% Extract the gamma data for the PR-6XX.
%
% This gets the measured spectrum and dark subtracts.
% It then computes the range of wavelengths around the peak
% that will be used when we compute power for each gamma band
% and computes the fraction of max power using regression on
% the spectra from the selected wavelength bands.
%
% Note: raw.gamma.rad is a structure array.  Each element
% of the structure array contains the measurements for one
% of the measured gamma bands.
if (size(cal.raw.gamma.rad,2) ~= cal.describe.nGammaBands)
    error('Mismatch between specified number of gamma bands and size of measurement struct array');
end
for k = 1:cal.describe.nGammaBands
    % Dark subtract with time correction.  As with primary spectra above,
    % there are two distinct ways we can do this, depending on whether we
    % calibrated around dark or around a specified background.  The end
    % result of this makes the two sets of measurements equivalent going
    % onwards.
    gammaTemp = cal.raw.gamma.rad(k).meas;
    gammaMeas{k} = gammaTemp - cal.computed.pr650MeanDark(:,ones(1,size(cal.raw.gamma.rad(k).meas,2)));
    [~,peakWlIndices{k}] = max(gammaMeas{k}(:,end));
    
    % Get wavelength indices around peak to use.
    cal.describe.peakWlIndex(k) = peakWlIndices{k}(1);
    cal.describe.minWlIndex(k) = cal.describe.peakWlIndex(k)-cal.describe.gammaNumberWlUseIndices;
    if (cal.describe.minWlIndex(k) < 1)
        cal.describe.minWlIndex(k) = 1;
    end
    cal.describe.maxWlIndex(k) = cal.describe.peakWlIndex(k)+cal.describe.gammaNumberWlUseIndices;
    if (cal.describe.maxWlIndex(k) > cal.describe.S(3))
        cal.describe.maxWlIndex(k) = cal.describe.S(3);
    end
    
    % Little check and then get power for each measured level for this measured band
    if (size(cal.raw.gamma.rad(k).meas,2) ~= cal.describe.nGammaLevels)
        error('Mismatch between specified number of gamma levels and size of measurement array');
    end
    
    for i = 1:cal.describe.nGammaLevels
        wavelengthIndices = cal.describe.minWlIndex(k):cal.describe.maxWlIndex(k);
        cal.computed.gammaData1{k}(i) = gammaMeas{k}(wavelengthIndices,end)\ ...
            gammaMeas{k}(wavelengthIndices,i); %#ok<*AGROW>
        cal.computed.gammaRatios(k,i+1).wavelengths = cal.computed.commonWls(wavelengthIndices);
        cal.computed.gammaRatios(k,i+1).ratios = gammaMeas{k}(wavelengthIndices,i) ./ gammaMeas{k}(wavelengthIndices,end);
    end
    cal.computed.gammaRatios(k,1).wavelengths = cal.computed.gammaRatios(k,2).wavelengths;
    cal.computed.gammaRatios(k,1).ratios = 0*cal.computed.gammaRatios(k,2).ratios;
end

% Fit each indivdually measured gamma function to finely spaced real valued
% primary levels.
cal.computed.gammaInput = linspace(0,1,cal.describe.nGammaFitLevels)';
for k = 1:cal.describe.nGammaBands
    cal.computed.gammaTableMeasuredBands(:,k) = [0 ; cal.computed.gammaData1{k}'];
    cal.computed.gammaTableMeasuredBandsFit(:,k) = OLFitGamma(cal.computed.gammaInputRaw,cal.computed.gammaTableMeasuredBands(:,k),cal.computed.gammaInput,cal.describe.gammaFitType);
end

% Interpolate the measured bands out across all of the bands
for l = 1:cal.describe.nGammaFitLevels
    cal.computed.gammaTable(l,:) = interp1(cal.describe.gamma.gammaBands',cal.computed.gammaTableMeasuredBandsFit(l,:)',(1:cal.describe.numWavelengthBands)','linear','extrap')';
end

% Make each band's gamma curve monotonic
for b = 1:cal.describe.numWavelengthBands
    cal.computed.gammaTable(:,b) = MakeMonotonic(cal.computed.gammaTable(:,b));
end

% Average gamma measurements over bands
switch (cal.describe.whichAverageGamma)
    case 'median'
        cal.computed.gammaTableAvg = median(cal.computed.gammaTableMeasuredBandsFit,2);
    case 'middle'
        nMeasuredBands = size(cal.computed.gammaTableMeasuredBandsFit,2);
        middleIndex = round(nMeasuredBands/2);
        cal.computed.gammaTableAvg = cal.computed.gammaTableMeasuredBandsFit(:,middleIndex);
    otherwise
        error('Unknown value for whichAverage gamma provided.')
end
cal.computed.gammaTableAvg = MakeMonotonic(cal.computed.gammaTableAvg);



