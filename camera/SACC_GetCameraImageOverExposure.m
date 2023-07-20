% SACC_GetCameraImageOverExposure.
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

% Choose which data to load either '0714' or '0719'.
measureDate = '0719';

if strcmp(measureDate,'0719')
    numTrial = '1st';
end

%% Load all images here.
nImages = length(exposureOptions);
for ee = 1:nImages
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration','Exposure',measureDate);
        switch measureDate
            case '0714'
                testFilename = GetMostRecentFileName(testFiledir, append(num2str(targetCyclePerDeg),'cpd_crop_',num2str(exposureOptions(ee))));
            case '0719'
                testFiledir = fullfile(testFiledir,numTrial);
                testFilename = GetMostRecentFileName(testFiledir,append(num2str(exposureOptions(ee)),'mm_crop'));
        end
        % We save all images here. The array looks like {dataType,
        % channel, SF}.
        images{ee} = imread(testFilename);
    else
        error('Cannot find data file');
    end
end

%% Plot the camera images.
PLOTIMAGE = false;
if (PLOTIMAGE)
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    for ee = 1:nImages
        subplot(3,3, ee);
        imshow(images{ee});
        title(sprintf('Exposure = %s mm',num2str(exposureOptions(ee))),'FontSize',15);
    end
end

%% Plot the sliced images.
PLOTSLICEDIMAGE = true;
if strcmp(measureDate,'0714')
    if (PLOTSLICEDIMAGE)
        
        % Make a new figure per each channel.
        figure;
        figurePosition = [0 0 1500 800];
        set(gcf,'position',figurePosition);
        
        % Add a grand title of the figure.
        sgtitle('Image contrast vs. Camera Exposure','FontSize',15);
        
        % Set the minimum distance between two peaks.
        minPeakDistance = 4;
        
        % Calculate the contrasts here.
        %
        % We will skip the first image with exposure = 10 mm, which is fully
        % closed, where we cannot measure the contrast.
        for ee = 1 : nImages
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
end

%% Plot constrasts over exposure level.
%
% Make a new figure.
if strcmp(measureDate,'0714')
    figure; hold on;
    xticks(exposureOptions);
    xlabel('Exposure (mm)','fontsize',15);
    ylabel('Mean Contrast','fontsize',15);
    title('Image contrasts over camera exposure','fontsize',15)
    ylim([0.5 1.05]);
    plot(exposureOptions,cell2mat(meanContrasts),'r.','markersize',11);
    errorbar(exposureOptions,cell2mat(meanContrasts),cell2mat(stdErrorContrasts),'r');
    e = get(gca,'children');
    e(1).LineStyle = 'none';
end

%% Plot the mean power over exposure level.
%
% Get mean power of the images.
for ii = 1:nImages
    meanPowerImages(ii) = mean(images{ii},'all');
end
figure; hold on;
xticks(exposureOptions);
xlabel('Exposure (mm)','fontsize',15);
ylabel('Mean power (dRGB)','fontsize',15);
plot(exposureOptions,meanPowerImages,'ro','markersize',11,'markerfacecolor','r','markeredgecolor','k');
