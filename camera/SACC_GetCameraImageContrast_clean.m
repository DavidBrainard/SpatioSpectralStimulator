% SACC_GetImageContrastWithFunction.
%
% It calculates the image contrasts using the function.
%
% See also:
%    SACC_GetImageContrastWithFunction

% History:
%    06/13/23   smo    - Wrote it.
%    06/27/23   smo    - Now we load all images and analyze it together.
%    07/13/23   smo    - Added a plot to compare the actual peaks of the
%                        spectrums used in each projector.
%    08/14/23   smo    - Made a clean version.

%% Initialize.
clear; close all;

%% Set variables.
%
% Initial measurements were made on 0613.
targetCyclePerDeg = {18};
projectorSettings = {'Raw'};
measureDate = '0814';
measureDistanceOptions = {'1d','1.66d','2d'};

% %% Plot the spectrum used.
% if (ispref('SpatioSpectralStimulator','SACCMaterials'))
%     testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
%     testFiledir = fullfile(testFiledir,'Calibration');
%     testFilename = GetMostRecentFileName(testFiledir,'SACCPrimary1');
%
%     % We save all images here. The array looks like {dataType,
%     % channel, SF}.
%     data = load(testFilename);
% end
%
% % Here we read out both old and new projector calibration files. For new
% % projector, we read the most recent one. For old projector, we read the
% % last one measured which is stored in 17th of the calibration file in
% % SACCPriamry1.
% idxFileOldProjector = 17;
% calData_oldProjector = data.cals{idxFileOldProjector};
% calData_newProjector = data.cals{end};
% S = calData_newProjector.rawData.S;
% wls = SToWls(S);
% spds_newProjector = calData_newProjector.processedData.P_device;
% spds_oldProjector = calData_oldProjector.processedData.P_device;
%
% % Get spds of the channels used.
% %
% % We used the same channels for both old and new projector [2 5 9 12 15] on
% % the date of 0613. We will set these numbers differently so that the
% % actual peaks are matched each other.
% switch measureDate
%     case '0613'
%         numChannelUsed_newProjector = [2 5 9 12 15];
%         numChannelUsed_oldProjector = [2 5 9 12 15];
%     case '0714'
%         numChannelUsed_newProjector = [1 3 7 12 15];
%         numChannelUsed_oldProjector = [1 3 7 11 14];
%     case '0718'
%         numChannelUsed_newProjector = [1 3 7 10 15];
%         numChannelUsed_oldProjector = [1 3 7 9 14];
%     case '0719'
%         numChannelUsed_newProjector = [1 3 7 10 15];
%         numChannelUsed_oldProjector = [1 3 7 9 14];
%     case '0814'
%         numChannelUsed_newProjector = [7];
%         numChannelUsed_oldProjector = [7];
% end
%
% spdsUsed_newProjector = spds_newProjector(:,numChannelUsed_newProjector);
% spdsUsed_oldProjector = spds_oldProjector(:,numChannelUsed_oldProjector);
%
% % Peaks of each spectrums.
% peaksUsed_newProjector = FindPeakSpds(spdsUsed_newProjector,'verbose',false);
% peaksUsed_oldProjector = FindPeakSpds(spdsUsed_oldProjector,'verbose',false);
%
% % Plot it.
% figure;
% figPosition = [0 0 1500 400];
% set(gcf,'position',figPosition);
% sgtitle('Subprimary channels used for measuring chromatic aberration','fontsize',15);
%
% % New projector.
% subplot(1,3,1); hold on;
% plot(wls,spds_newProjector,'k-','linewidth',0.8);
% plot(wls,spdsUsed_newProjector,'-','color',[1 0 0 0.3],'linewidth',4);
% xlabel('Wavelength (nm)','fontsize',15);
% ylabel('Spectral power','fontsize',15);
% xticks([380:80:780]);
% xticklabels([380:80:780]);
% ylim([0 max(max(spds_newProjector))*1.05])
% f = get(gca, 'children');
% legend(f(flip([1 17])),'All channels','Used channels',...
%     'location','northwest','fontsize',13)
% title('SACCSFA','fontsize',15);
% subtitle(sprintf('Peaks = (%s) nm',num2str(peaksUsed_newProjector)),'fontsize',14);
%
% % Old projector.
% subplot(1,3,2); hold on;
% plot(wls,spds_oldProjector,'k-','linewidth',0.8);
% plot(wls,spdsUsed_oldProjector,'-','color',[0 1 0 0.3],'linewidth',4);
% xlabel('Wavelength (nm)','fontsize',15);
% ylabel('Spectral power','fontsize',15);
% xticks([380:80:780]);
% xticklabels([380:80:780]);
% ylim([0 max(max(spds_newProjector))*1.05])
% f = get(gca, 'children');
% legend(f(flip([1 17])),'All channels','Used channels',...
%     'location','northwest','fontsize',13)
% title('Raw (Old Projector)','fontsize',15);
% subtitle(sprintf('Peaks = (%s) nm',num2str(peaksUsed_oldProjector)),'fontsize',14);
%
% % Comparison New vs. Old projector.
% subplot(1,3,3), hold on;
% plot(wls,spdsUsed_newProjector,'-','color',[1 0 0 0.3],'linewidth',4);
% plot(wls,spdsUsed_oldProjector,'-','color',[0 1 0 0.3],'linewidth',4);
% xlabel('Wavelength (nm)','fontsize',15);
% ylabel('Spectral power','fontsize',15);
% xticks([380:80:780]);
% xticklabels([380:80:780]);
% ylim([0 max(max(spds_newProjector))*1.05])
% f = get(gca, 'children');
% legend(flip(f([1 6])),'SACCSFA','Raw (Old Projector)',...
%     'location','northwest','fontsize',13)
% title('SACCSFA vs. Raw','fontsize',15);
%
% % Save the peak wavelengths in string. We will use this for legend in the
% % following plot.
% for pp = 1:length(peaksUsed_newProjector)
%     peakWls{pp} = num2str(peaksUsed_newProjector(pp));
% end

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
