% SACC_ContrastPrintedImage.
%
% It calculates the image contrasts using the function.
%
% See also:
%    SACC_GetCameraImageContrast

% History:
%    08/15/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
targetCyclePerDeg = {3,6,9,12,18};
measureDate = '0829';

%% Load all images here.
for cc = 1:8
    targetChOptions{cc} = append('Ch',num2str(cc));
end
nTargetChs = length(targetChOptions;

% Make a loop here.
for 
targetChTemp = targetChOptions
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDate,targetChTemp);
else
    error('Cannot find data file list!');
end

nSFs = length(targetCyclePerDeg);
for dd = 1:nSFs
    % Get the file name of the images.
    testFilenameTemp = GetMostRecentFileName(testFiledir,...
        append(num2str(targetCyclePerDeg{dd}),'cpd_crop'));
    
    % We save all images here.
    images{dd} = imread(testFilenameTemp);
end

%% Plot the camera images.
PLOTIMAGE = true;
if (PLOTIMAGE)
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    for dd = 1:nSFs
        subplot(nSFs,1,dd);
        imshow(images{dd});
        title(sprintf('%d cpd',targetCyclePerDeg{dd}))
    end
end

%% Plot the sliced images.
PLOTSLICEDIMAGE = true;
if (PLOTSLICEDIMAGE)
    % Make a new figure.
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    % Contrast calculation happens here.
    for dd = 1:nSFs
        % Set min distance between adjacent peaks.
        if dd == 1
            minPeakDistance = 35;
        elseif dd == 2
            minPeakDistance = 20;
        else
            minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        subplot(nSFs,1,dd);
        title(sprintf('%d cpd',targetCyclePerDeg{dd}),'fontsize',15);
        
        % Calculate contrasts.
        contrastsRawTemp = GetImgContrast(images{dd},'minPeakDistance',minPeakDistance);
        contrastsRaw{dd} = contrastsRawTemp;
        meanContrasts{dd} = mean(contrastsRawTemp);
        stdErrorContrasts{dd} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
    end
    
    % Save the plot if you want.
    SAVETHEPLOT = false;
    if (SAVETHEPLOT)
        testFileFormat = '.tiff';
        testFilenameTemp = sprintf('%s_%s',projectorSettings{dd},channels{cc});
        saveas(gcf,append(testFilenameTemp,testFileFormat));
        disp('Plot has been saved successfully!');
    end
end

%% Plot contrasts over spatial frequency.
figure; clf;
plot(cell2mat(targetCyclePerDeg),cell2mat(meanContrasts),...
    'ko-','markeredgecolor','k','markerfacecolor','b','markersize',10);
ylim([0 1]);
xlabel('Spatial Frequency (cpd)','fontsize',15);
ylabel('Mean Contrasts','fontsize',15);
xticks(cell2mat(targetCyclePerDeg));

%% Plot spectra.
if (strcmp(measureDate,'0829'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDate);
    testFilename = 'spd_combiLED.mat';
    spdData = load(fullfile(testFiledir,testFilename));
    
    % Extract the spd data and flip left to right, becasue the measurement was
    % done from Ch8 (high, 652 nm) to Ch1 (low, 406 nm).
    spd = spdData.spd;
    spd = fliplr(spd);
    S = [380 2 201];
    wls = SToWls(S);
    peaks_spd = FindPeakSpds(spd,'verbose',false);
    
    % Plot it - Raw.
    figure; clf;
    plot(wls,spd,'linewidth',1);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    for ll = 1:size(spd,2)
        legendHandles{ll} = append(num2str(peaks_spd(ll)),' nm');
    end
    legend(legendHandles,'fontsize',15);
    
    % Plot it - Normalized.
    figure; clf;
    plot(wls,spd./max(spd),'linewidth',1);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    legend(legendHandles,'fontsize',15);
end
