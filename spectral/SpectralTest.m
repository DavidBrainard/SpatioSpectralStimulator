% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

% Load in Excel file ad get raw individual LED spectra
projectName = 'SpatioSpectralStimulator';
projectRoot = tbLocateProject(projectName);
dataDir = 'data';
ledFilename = 'LUXEON_CZ_SpectralPowerDistribution.xlsx';
ledXlsData = xlsread(fullfile(projectRoot,dataDir,ledFilename));
nLeds = size(ledXlsData,2)/2;
index = 1;
for ii = 1:nLeds
    ledSpectraRaw{ii}(:,1) = ledXlsData(:,index);
    ledSpectraRaw{ii}(:,2) = ledXlsData(:,index+1);
    index = index+ 2;
end

%% Get rid of NaN's, sort, and fit each spectrum
%
% Then interpolate to standard wavelengh support to
% form basis matrix
smoothingParam = 0.5;
S = [380 1 401];
wls = SToWls(S);
for ii = 1:nLeds
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
    ledBasis(:,ii) = interp1(ledSmooth{ii}(:,1),ledSmooth{ii}(:,2),wls,'linear',0);
end

%% Plot basis
figure; clf; hold on
for ii = 1:nLeds
    theColor = rand(1,3);
    plot(ledSpectra{ii}(:,1),ledSpectra{ii}(:,2),'o','Color',theColor, ...
        'MarkerFaceColor',theColor,'MarkerSize',2);
    plot(wls,ledBasis(:,ii),'Color',theColor,'LineWidth',2);
end

%% Fit a spectrum with the basis
theTemp = 6500;
theSpd = GenerateBlackBody(theTemp,wls);
weights = lsqnonneg(ledBasis,theSpd);
theFitSpd = ledBasis*weights;
figure; clf; hold on
plot(wls,theSpd,'k','LineWidth',3);
plot(wls,theFitSpd,'r','LineWidth',2);

