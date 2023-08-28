% SACC_PlotSpectralCharacter.
%
% This gives some plots of SACCSFA spectral characterization.

% History:
%    08/27/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Read out the data.
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Calibration');
    testFilename = 'SACCPrimary1.mat';
    data = load(fullfile(testFiledir,testFilename));
else
    error('Cannot find data file list!');
end

% We will use the most recent measurement data.
data = data.cals{end};

%% Extract the data that we will use.
S = data.rawData.S;
wls = SToWls(S);
spd_primariesRaw = data.processedData.P_device;
spd_primariesNorm = spd_primariesRaw./max(spd_primariesRaw);

% Find peaks of each primary and sort in ascending order.
peaks_primaries = FindPeakSpds(spd_primariesRaw,'verbose',false);
[peaks_primaries_sorted i] = sort(peaks_primaries);

%% Plot it.
%
% Primary spectral - Raw.
figure; clf;
plot(wls,spd_primariesRaw,'linewidth',1);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spectral Power','fontsize',15);
ylim([0 max(spd_primariesRaw,[],'all')*1.1]);
xticks([380:80:780]);

% Primary spectral - Normalized.
figure; clf;
plot(wls,spd_primariesNorm,'linewidth',1);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spectral Power','fontsize',15);
ylim([0 1]);
xticks([380:80:780]);
f = get(gca,'children');
f = flip(f);
for pp = 1:length(f)
    lineColorOptions(pp,:) = f(pp).Color;
end

%% Gamma function.
figure; clf; hold on;
nPrimaries = size(spd_primariesRaw,2);
gammaInput = data.rawData.gammaInput; 
for pp = 1:nPrimaries
    lineColorTemp = lineColorOptions(pp,:);
    gammaTableTemp = data.rawData.gammaTable(:,pp)';
    
    % Plot it here.
    plot([0 gammaInput],[0 gammaTableTemp],'o-',...
        'markerfacecolor',lineColorTemp,'color',lineColorTemp);
    
    % Linear fit per each subprimary to get the slope.
    pFit(:,pp) = polyfit(gammaInput,gammaTableTemp,1); 
end
xlabel('Settings input','fontsize',15);
ylabel('Settings output','fontsize',15);
ylim([0 1]);

%% CIE xy coordinates.
%
% Load color matching functions.
load T_xyzJuddVos
T_XYZ = T_xyzJuddVos';
T_XYZ = SplineCmf(S_xyzJuddVos,T_XYZ',S,2);

% Get xy coordinates.
XYZ = T_XYZ*spd_primariesRaw;
xyY = XYZToxyY(XYZ);

% Spectral locus.
spectralLocus = XYZToxyY(T_XYZ);
spectralLocus(:,end+1) = spectralLocus(:,1);

% Standard color space.
xyz_sRGB = [0.64 0.30 0.15; 0.33 0.60 0.06; 0.2126 0.7152 0.0722];
xyz_sRGB(:,end+1) = xyz_sRGB(:,1);

% Plot it.
figure; clf; hold on;

% sRGB space.
plot(xyz_sRGB(1,:),xyz_sRGB(2,:),'-','linewidth',2,'color',[0 0 0 0.3]);

% Connect the points between subprimaries to draw the boundary.
xyY_sorted = xyY(:,i);
xyY_sorted(:,end+1) = xyY_sorted(:,1);
plot(xyY_sorted(1,:),xyY_sorted(2,:),'-','color',[0 0 1 0.3],'linewidth',2);

% Plot the data points.
for pp = 1:nPrimaries
    plot(xyY(1,pp),xyY(2,pp),'o','markerfacecolor',lineColorOptions(pp,:),'markeredgecolor','k','markersize',8);
end

% Calculate the area of sRGB and SACCSFA gamut.
area_sRGB = polyarea(xyz_sRGB(1,:), xyz_sRGB(2,:));
area_SACCSFA = polyarea(xyY_sorted(1,:), xyY_sorted(2,:));
percent_SACCSFA2sRGB = area_SACCSFA/area_sRGB;

% Spectral locus.
plot(spectralLocus(1,:),spectralLocus(2,:),'k-');
xlabel('CIE x','fontsize',15);
ylabel('CIE y','fontsize',15);
legend('sRGB','SACCSFA','fontsize',15);

%% Spectral shifts.
%
% Read out the gamma curve spectra.
% Array is sorted as channel (16) x gamma input level (10) x number of points (201).
spd_gammaCurve = data.rawData.gammaCurveMeanMeasurements;
nGammaInputs = length(data.rawData.gammaInput);

% Plot it.
figure; clf; hold on;
figurePosition = [0 0 1000 1000];
set(gcf, 'position', figurePosition);

% Sort the spectra and line color options as well.
spd_gammaCurve_sorted = spd_gammaCurve(i,:,:);
lineColorOptions_sorted = lineColorOptions(i,:);

% Find a max value for setting the y-axis of the following plot.
max_ylim = max(imresize(spd_gammaCurve,[size(spd_gammaCurve,1),size(spd_gammaCurve,2)*size(spd_gammaCurve,3)]),[],'all');

% Make a loop for the spectra plot.
for pp = 1:nPrimaries
    subplot(4,4,pp); hold on;
    title(sprintf('%d nm',peaks_primaries_sorted(pp)),'fontsize',15);
    for gg = 1:nGammaInputs
        spd_temp = squeeze(spd_gammaCurve_sorted(pp,gg,:));
        plot(wls,spd_temp,'color',lineColorOptions_sorted(pp,:));
        xlabel('Wavelength (nm)');
        ylabel('Spectral power');
        xticks([380:80:780]);
        ylim([0 max_ylim]);
    end
end
