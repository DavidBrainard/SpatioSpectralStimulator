% SpectraOverFilters.
%
% This shows the spectra of test images over using different filters, A, B,
% C, D, and E.

% History:
%    6/30/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
filterOptions = {'A', 'B', 'C', 'D', 'E'};

%% Load the test image spectra.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCData');
    testFiledir = fullfile(testFiledir,'CheckCalibration');
    testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck');
    spdTestData = load(testFilename);
else
    error('Cannot find the data file!');
end

% Get spectrum range.
S = spdTestData.theData.S;
wls = SToWls(S);

% Read out the spectra.
spd_BG = spdTestData.screenBgSpd;
spd_Test = spdTestData.ptCldScreenSpdMeasuredCheckCal;

% Image contrasts.
testContrasts = spdTestData.theData.rawMonochromeUnquantizedContrastCheckCal .* spdTestData.spatialGaborTargetContrast;

% Get red and green spectra on the test image. Here we will choose the
% highest contrast possible within the range.
idxMaxContrast = find(testContrasts == max(testContrasts));
idxMinContrast = find(testContrasts == min(testContrasts));
spd_Red = spd_Test(:,idxMaxContrast);
spd_Green = spd_Test(:,idxMinContrast);

% Plot it.
figure; hold on;
plot(wls,spd_BG,'k-','linewidth',2);
plot(wls,spd_Red,'r-','linewidth',2);
plot(wls,spd_Green,'g-','color',[0.1 0.7 0.1],'linewidth',2);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spectral power','fontsize',15);
legend('Background', sprintf('Red (Image Contrast = %.2f)', testContrasts(idxMaxContrast)),...
    sprintf('Green (Image Contrast = %.2f)', testContrasts(idxMinContrast)),'fontsize',12);

%% Load the filter data.
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Filters');
    testFilename = fullfile(testFiledir,'GoggleSpectraSummary_UPenn.xlsx');
    spdFilters = readtable(testFilename,'sheet','SPD');
else
    error('Cannot find the data file!');
end

% Get wavelength range.
S_filter = spdFilters.Wave;

% Get SPDs of the filters.
nFilters = length(filterOptions);
for ff = 1:nFilters
    spds(:,ff) = table2array(spdFilters(:,ff+1));
end

% Match the wavelength.
spd_Filters = SplineRaw(S_filter,spds,S);

% Plot it.
figure; hold on;
for ff = 1:nFilters
    plot(wls,spd_Filters(:,ff), 'color', [0.2*ff 0.2*ff 0], 'linewidth', 2);
end
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Transmittance','fontsize',15);
legend(filterOptions,'location','southeast','fontsize',13);

%% Calculate the spectra (test image x filters).
%
% New figure.
fig=figure; hold on;
figSize = 700;
figPosition = [0 0 figSize figSize];
set(fig, 'position', figPosition);

% Set test spectra.
testSpectra = {spd_BG, spd_Red, spd_Green};
testSpectraName = {'Background', 'Red', 'Green'};

% Background over the filters.
for cc = 1:length(testSpectra)
    % Make a subplot per each test spectra.
    subplot(3,1,cc); hold on;
    
    % Raw spectra without using filters.
    plot(wls, testSpectra{cc}.*max(max(spd_Filters)), 'k-','color', [0 0 0 0.4],'linewidth', 6);
    
    % Calculate the spectra over the filters here.
    for ff = 1:nFilters
        plot(wls, testSpectra{cc}.*spd_Filters(:,ff),...
            'color', [0.2*ff 0.2*ff 0], 'linewidth', 2);
    end
    xlabel('Wavelength (nm)');
    ylabel('Spectral Power');
    ylim([0 0.3]);
    legend({'Raw' filterOptions{:}});
    title(testSpectraName{cc});
end

%% Calculate the cone contrasts.
%
% Load cone sensitivities.
T_cones = spdTestData.theData.T_cones;

% Get all spectra with filters.
for ff = 1:nFilters
    spd_BG_AllFilters(:,ff) = spd_BG .* spd_Filters(:,ff);
    spd_Red_AllFilters(:,ff) = spd_Red .* spd_Filters(:,ff);
    spd_Green_AllFilters(:,ff) = spd_Green .* spd_Filters(:,ff);
end

% Cone excitations.
exicitations_BG_Raw = T_cones * spd_BG;
excitations_Red_Raw = T_cones * spd_Red;
excitations_Green_Raw = T_cones * spd_Green;

excitations_BG_Filters = T_cones * spd_BG_AllFilters;
excitations_Red_Filters = T_cones * spd_Red_AllFilters;
excitations_Green_Filters = T_cones * spd_Green_AllFilters;

% Calculate contrasts.
contrastsRaw_Red = ExcitationToContrast(excitations_Red_Raw,exicitations_BG_Raw);
contrastsRaw_Green = ExcitationToContrast(excitations_Green_Raw,exicitations_BG_Raw);
for ff = 1:nFilters
    contrastsFilters_Red(:,ff) = ExcitationToContrast(excitations_Red_Filters(:,ff), excitations_BG_Filters(:,ff));
    contrastsFilters_Green(:,ff) = ExcitationToContrast(excitations_Green_Filters(:,ff), excitations_BG_Filters(:,ff));
end

% Merge all contrast arrays with and without filters.
contrastsAll_Red = [contrastsRaw_Red contrastsFilters_Red];
contrastsAll_Green = [contrastsRaw_Green contrastsFilters_Green];

% Plot it.
coneOptions = {'L', 'M', 'S'};
barColor = {'k', [0.1 0.7 0.1], 'b'};
xaxisHandles = {'Raw', filterOptions{:}};

% RED.
figure; hold on;
figPosition = [0 0 figSize+300 figSize];
set(gcf,'position',figPosition);

for cc = 1:length(coneOptions)
    subplot(1,3,cc); hold on;
    
    % By filter.
    for ff = 1:nFilters+1
        bar(ff, contrastsAll_Red(cc,ff), 'FaceColor',[0.2*(ff-1) 0.2*(ff-1) 0],'EdgeColor',[0 0 0]);
    end
    ylim([-0.3 0.25]);
    xticks([1:1:6]);
    xticklabels(xaxisHandles)
    xlabel('Filters','fontsize',14);
    ylabel('Cone Contrast','fontsize',14);
    title(append(coneOptions{cc},' cone'),'fontsize',14);
end
sgtitle('Cone contrast changes over the filters - (RED)','fontsize',15);

% Green.
figure; hold on;
figPosition = [0 0 figSize+300 figSize];
set(gcf,'position',figPosition);

for cc = 1:length(coneOptions)
    subplot(1,3,cc); hold on;
    
    % By filter.
    for ff = 1:nFilters+1
        bar(ff, contrastsAll_Green(cc,ff), 'FaceColor',[0.2*(ff-1) 0.2*(ff-1) 0],'EdgeColor',[0 0 0]);
    end
    ylim([-0.3 0.25]);
    xticks([1:1:6]);
    xticklabels(xaxisHandles)
    xlabel('Filters','fontsize',15);
    ylabel('Cone Contrast','fontsize',15);
    title(append(coneOptions{cc},' cone'),'fontsize',15);
end
sgtitle('Cone contrast changes over the filters - (Green)','fontsize',15);
