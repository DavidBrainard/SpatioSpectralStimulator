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
xlabel('Wavelength (nm)');
ylabel('Spectral power');
legend('Background', sprintf('Red (Image Contrast %.2f)', testContrasts(idxMaxContrast)),...
    sprintf('Green (Image Contrast %.2f)', testContrasts(idxMinContrast)));

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
xlabel('Wavelength (nm)');
ylabel('Transmittance');
legend(filterOptions);

%% Calculate the spectra (test image x filters).
%
% New figure.
fig=figure; hold on;
figSize = 700;
set(fig, 'position', [0 0 figSize figSize]);

% Set test spectra.
testSpectra = {spd_BG, spd_Red, spd_Green};
testSpectraName = {'Background', 'Red', 'Green'};

% Background over the filters.
for tt = 1:length(testSpectra)
    % Make a subplot per each test spectra.    
    subplot(3,1,tt); hold on;
    
    % Raw spectra without using filters.
    plot(wls, testSpectra{tt}.*max(max(spd_Filters)), 'k-','color', [0 0 0 0.4],'linewidth', 6);
    
    % Calculate the spectra over the filters here.
    for ff = 1:nFilters
        plot(wls, testSpectra{tt}.*spd_Filters(:,ff),...
            'color', [0.2*ff 0.2*ff 0], 'linewidth', 2);
    end
    xlabel('Wavelength (nm)');
    ylabel('Spectral Power');
    ylim([0 0.3]);
    legend({'Raw' filterOptions{:}});
    title(testSpectraName{tt});
end
