function theBasis = sssGetSpectralBasis(S,varargin)
% Get a simulated spectral basis for the spatial-spectral simulator
%
% Syntax:
%   theBasis = sssGetSpectralBasis(S)
%
% Description:
%    Get a simulated spectral basis for the spatial-spectral simulator.
%    This is useful for testing and design development.
%
%    Default values based on digitized LED spectra plus some narrowband
%    Gaussians.
%
% Inputs:
%    S          - Wavelength support in PTB format.  May be S or wls
%
% Outputs:
%    theBasis   - Basis matrix, with the simulated spectral basis in the
%                 columns.
%
% Optional key/value pairs:
%    plotBasis  - Boolean (default false).  Make a plot of the basis.
%
% See Also:
%   
%

% History:
%    05/03/20  dhb  Pulled out from test script.   


% Examples:
%{
    S = [380 1 401];
    B = sssGetSpectralBasis(S,'plotBasis',true);
%}

% Massage varargin and parse
p = inputParser;
p.addParameter('plotBasis', false, @islogical);
p.parse(varargin{:});

% Put S in desired format
S = MakeItS(S);

% Load in Excel file ad get raw individual LED spectra
projectName = 'SpatioSpectralStimulator';
projectRoot = tbLocateProject(projectName,[],'verbose',false);
dataDir = 'data';
ledFilename = 'LUXEON_CZ_SpectralPowerDistribution.xlsx';
ledXlsData = xlsread(fullfile(projectRoot,dataDir,ledFilename));
nLedsData = size(ledXlsData,2)/2;
index = 1;
for ii = 1:nLedsData
    ledSpectraRaw{ii}(:,1) = ledXlsData(:,index);
    ledSpectraRaw{ii}(:,2) = ledXlsData(:,index+1);
    index = index+ 2;
end

%% Get rid of NaN's, sort, and fit each spectrum
%
% Then interpolate to standard wavelengh support to
% form basis matrix
smoothingParam = 0.5;
wls = SToWls(S);
for ii = 1:nLedsData
    index = ~isnan(ledSpectraRaw{ii}(:,1));
    ledSpectra{ii}(:,1) = ledSpectraRaw{ii}(index,1);
    ledSpectra{ii}(:,2) = ledSpectraRaw{ii}(index,2);  
    [~,index] = sort(ledSpectra{ii}(:,1));
    ledSpectra{ii}(:,1) = ledSpectra{ii}(index,1);
    ledSpectra{ii}(:,2) = ledSpectra{ii}(index,2);
    
    % Fit with smooth spline
    fitobj = fit(ledSpectra{ii}(:,1),ledSpectra{ii}(:,2),'smoothingspline', ...
        'SmoothingParam',smoothingParam);
    ledSmooth{ii}(:,1) = linspace(ledSpectra{ii}(1,1),ledSpectra{ii}(end,1),1000);
    ledSmooth{ii}(:,2) = feval(fitobj,ledSmooth{ii}(:,1));
    ledSmooth{ii}(:,2) = ledSmooth{ii}(:,2)/max(ledSmooth{ii}(:,2));
    ledBasisData(:,ii) = interp1(ledSmooth{ii}(:,1),ledSmooth{ii}(:,2),wls,'linear',0);
end

%% Take out broadband
whichToRemove = [6 7 8];
whichToKeep = setdiff(1:nLedsData,whichToRemove);
ledBasis = ledBasisData(:,whichToKeep);
nLeds = size(ledBasisData,2);

%% Add monochromatic
fullWidthHalfMax = 15;
gaussVar = FWHMToStd(fullWidthHalfMax)^2;
monochromaticToAdd = [440 540];
gaussBasis = MakeGaussBasis(wls,monochromaticToAdd,gaussVar*ones(size(monochromaticToAdd)));
gaussBasis = gaussBasis/max(gaussBasis(:));
theBasis = [ledBasis, gaussBasis];
nBasis = size(theBasis,2);

%% Plot basis
if (p.Results.plotBasis)
    figure; clf; hold on
    for ii = 1:nLedsData
        theColor = rand(1,3);
        plot(ledSpectra{ii}(:,1),ledSpectra{ii}(:,2),'o','Color',theColor, ...
            'MarkerFaceColor',theColor,'MarkerSize',2);
        plot(wls,ledBasisData(:,ii),'Color',theColor,'LineWidth',2);
    end
    
    figure; clf; hold on
    for ii = 1:nBasis
        theColor = rand(1,3);
        plot(wls,theBasis(:,ii),'Color',theColor,'LineWidth',2);
    end
end

