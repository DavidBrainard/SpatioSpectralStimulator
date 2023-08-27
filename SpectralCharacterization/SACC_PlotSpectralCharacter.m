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
for pp = 1:nPrimaries
    lineColorTemp = lineColorOptions(pp,:);
    plot([0 data.rawData.gammaInput],[0 data.rawData.gammaTable(:,pp)'],'ko-',...
        'markerfacecolor',lineColorTemp,'color',lineColorTemp);
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
plot(xyz_sRGB(1,:),xyz_sRGB(2,:),'b-','linewidth',1.5,'color',[0 0 1 0.3]);
for pp = 1:nPrimaries
    plot(xyY(1,pp),xyY(2,pp),'ro','markerfacecolor',lineColorOptions(pp,:),'markeredgecolor','k','markersize',8);
end 
plot(spectralLocus(1,:),spectralLocus(2,:),'k-');
xlabel('CIE x','fontsize',15);
ylabel('CIE y','fontsize',15);
legend('sRGB','fontsize',15);

%% Spectral shifts.
figure; clf; hold on;
