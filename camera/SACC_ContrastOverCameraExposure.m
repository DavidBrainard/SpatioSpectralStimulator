% SACC_ContrastOverCameraExposure.
%
% It calculates the image contrasts over different camera exposure based on
% the measured images.
%
% See also:
%    SACC_GetCameraImageContrast.

% History:
%    07/14/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
targetCyclePerDeg = 18;
exposureOptions = [10:2.5:30];
measureDate = '0714';

%% Load all images here.
nExposures = length(exposureOptions);
for ee = 1:nExposures
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration','Exposure',measureDate);
        testFilename = GetMostRecentFileName(testFiledir, append(num2str(targetCyclePerDeg),'cpd_crop_',num2str(exposureOptions(ee))));
        
        % We save all images here. The array looks like {dataType,
        % channel, SF}.
        images{ee} = imread(testFilename);
    else
        error('Cannot find data file');
    end
end

%% Plot the camera images.
PLOTIMAGE = true;
if (PLOTIMAGE)
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    for ee = 1:nExposures
        subplot(3,3, ee);
        imshow(images{ee});
        title(sprintf('Exposure = %s mm',num2str(exposureOptions(ee))),'FontSize',15);
    end
end

%% Plot the sliced images.
PLOTSLICEDIMAGE = true;
if (PLOTSLICEDIMAGE)
    
    % Make a new figure per each channel.
    figure;
    figurePosition = [0 0 1500 800];
    set(gcf,'position',figurePosition);
    
    % Add a grand title of the figure.
    sgtitle('Contrast vs. Camera Exposure','FontSize',15);
    
    % Set the minimum distance between two peaks.
    minPeakDistance = 4;
    
    % Calculate the contrasts here.
    %
    % We will skip the first image with exposure = 10 mm, which is fully
    % closed, where we cannot measure the contrast.
    for ee = 1 : nExposures
        if ee == 1
            contrastsRaw{ee} = 0;
            meanContrasts{ee} = 0;
            stdErrorContrasts{ee} = 0;
            continue
        end
        subplot(3,3,ee);
        contrastsRawTemp = GetImgContrast(images{ee},'minPeakDistance',minPeakDistance);
        contrastsRaw{ee} = contrastsRawTemp;
        meanContrasts{ee} = mean(contrastsRawTemp);
        stdErrorContrasts{ee} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
        title(sprintf('Exposure = %s mm',num2str(exposureOptions(ee))),'FontSize',15);
        ylim([0 16]);
    end
    
    % Save the plot if you want.
    SAVETHEPLOT = false;
    if (SAVETHEPLOT)
        testFileFormat = '.tiff';
        testFilename = sprintf('%s_%s',num2str(exposureOptions(ee)));
        saveas(gcf,append(testFilename,testFileFormat));
        disp('Plot has been saved successfully!');
    end
end

%% Plot constrasts over exposure level.
%
% Make a new figure.
figure; hold on;
xticks(exposureOptions);
xlabel('Exposure (mm)');
ylabel('Mean Contrast');
title('Contrasts over exposure')
ylim([0.5 1.05]);
plot(exposureOptions,cell2mat(meanContrasts),'r.','markersize',11);
errorbar(exposureOptions,cell2mat(meanContrasts),cell2mat(stdErrorContrasts),'r');
e = get(gca,'children');
e(1).LineStyle = 'none';
