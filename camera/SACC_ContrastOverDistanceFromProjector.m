% SACC_ContrastOverDistanceFromProjector.
%
% It calculates the image contrasts using the function.
%
% See also:
%    SACC_GetImageContrastWithFunction

% History:
%    08/14/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
%
% Initial measurements were made on 0613.
targetCyclePerDeg = {18};
projectorSettings = {'Raw'};
measureDate = '0814';
measureDistanceOptions = {'1d','1.66d','2d'};

%% Load all images here.
nSFs = length(targetCyclePerDeg);
nProjectorSettings = length(projectorSettings);

% Data type.
for dd = 1:length(projectorSettings)
    projectorSettingTemp = projectorSettings{dd};
    
    % Get channel name from the existing folders.
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        if dd == 3
            testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration',measureDate,projectorSettingTemp,'Focus Separately');
        else
            testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration',measureDate,projectorSettingTemp);
        end
        
        testFileList = dir(fullfile(testFiledir,'Ch*'));
    else
        error('Cannot find data file list!');
    end
    
    % Make a loop for Channel and Spatial frequency.
    nChannels = length(testFileList);
    for cc = 1:nChannels
        channels{cc} = testFileList(cc).name;
        
        % Extract only numbers. We are going to sort the array in an
        % ascending order.
        numChannelTemp = regexp(channels{cc}, '\d+', 'match');
        numChannels(cc) = str2double(numChannelTemp);
    end
    % Sorting the array (double array).
    [numChannelsSorted i] = sort(numChannels,'ascend');
    
    % We sort the channels in an ascending order (string array).
    channelsSorted = channels(i);
    
    for cc = 1:nChannels
        channelTemp = channelsSorted{cc};
        
        for dd = 1:nSFs
            % Get the file name of the images.
            testFiledirTemp = fullfile(testFiledir,channelTemp);
            testFilename = GetMostRecentFileName(testFiledirTemp,append(num2str(targetCyclePerDeg{dd}),'cpd_'));
            
            % We save all images here. The array looks like {dataType,
            % channel, SF}.
            images{dd,cc,dd} = imread(testFilename);
        end
    end
end

%% Load images temp.
%
% Get the file name of the images.
clear images;

nMeasureDistanceOptions = length(measureDistanceOptions);
for dd = 1:nMeasureDistanceOptions
    testFiledirTemp = fullfile(testFiledir,channelTemp);
    testFilename = GetMostRecentFileName(testFiledirTemp,...
        append(num2str(targetCyclePerDeg{:}),'cpd_',measureDistanceOptions{dd}));
    
    % We save all images here. The array looks like {dataType,
    % channel, SF}.
    images{dd} = imread(testFilename);
end

%% Plot the camera images.
PLOTIMAGE = true;
if (PLOTIMAGE)
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    for dd = 1:nMeasureDistanceOptions
        subplot(          nMeasureDistanceOptions,1,dd);
        imshow(images{dd});
        title(sprintf('%s',measureDistanceOptions{dd}))
    end
end

%% Plot the sliced images.
PLOTSLICEDIMAGE = true;
if (PLOTSLICEDIMAGE)
    for dd = 1:nProjectorSettings
        for cc = 1:nChannels
            % Make a new figure per each channel.
            figure;
            figurePosition = [0 0 800 800];
            set(gcf,'position',figurePosition);
            
            % Add a grand title of the figure.
            peakWls = {'530'};
            sgtitle(sprintf('%s - %s (%s nm)',projectorSettings{dd},channels{cc},peakWls{cc}),'FontSize',15);
            
            minPeakDistance = 5;
            % Calculate the contrasts here.
            for dd = 1:nMeasureDistanceOptions
                subplot(nMeasureDistanceOptions,1,dd);
                title(sprintf('%s',measureDistanceOptions{dd}),'fontsize',15);
                contrastsRawTemp = GetImgContrast(images{dd},'minPeakDistance',minPeakDistance);
                contrastsRaw{dd} = contrastsRawTemp;
                meanContrasts{dd} = mean(contrastsRawTemp);
                stdErrorContrasts{dd} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
            end
            
            % Save the plot if you want.
            SAVETHEPLOT = false;
            if (SAVETHEPLOT)
                testFileFormat = '.tiff';
                testFilename = sprintf('%s_%s',projectorSettings{dd},channels{cc});
                saveas(gcf,append(testFilename,testFileFormat));
                disp('Plot has been saved successfully!');
            end
        end
    end
end

%% Contrasts over the distance from the projector.
figure; clf;
numMeasureDistanceOptions = [1 1.66 2];
plot(numMeasureDistanceOptions,cell2mat(meanContrasts),'ko','markeredgecolor','k','markerfacecolor','b','markersize',10);
ylim([0 1]);
xlabel('Distance from the projector (no unit)','fontsize',15);
ylabel('Mean Contrasts','fontsize',15);
xticks(numMeasureDistanceOptions);
