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
measureDate = '0814';

%% Load all images here.
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDate);
else
    error('Cannot find data file list!');
end

nSFs = length(targetCyclePerDeg);
for dd = 1:nSFs
    % Get the file name of the images.
    testFilenameTemp = GetMostRecentFileName(testFiledir,...
        append(num2str(targetCyclePerDeg{dd}),'cpd_focused_crop'));
    
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
