% GetImageContrastWithFunction.
%
% It calculates the image contrasts using the function.

% History:
%    06/13/23   smo    - Wrote it.
%    06/27/23   smo    - Now we load all images and analyze it together.

%% Initialize.
clear; close all;

%% Set variables.
%
% For the data type, 'set1' is raw measurement, 'set2' is on the SACCSFA
% optical system.
%
% Channels corresponds, respectively, to 448, 506, 558, 618, and 658 nm
% peaks.
%
% All measurements were made on 0613.
targetCyclePerDeg = {3, 6, 9, 12, 18};
dataTypes = {'set1', 'set2'};
dataTypeSetting = {'Raw','SACCSFA'};
measureDate = '0613';
channels = {'Ch2', 'Ch5', 'Ch9', 'Ch12', 'Ch15'};
peakWls = {'448', '506', '558', '618', '658'};

%% Load all images here.
nChannels = length(channels);
nSFs = length(targetCyclePerDeg);
nDataTypes = length(dataTypes);

% Data type.
for dd = 1:length(dataTypes)
    dataType = dataTypes{dd};
    % Channel.
    for cc = 1:nChannels
        whichChannel = channels{cc};
        % Spatial frequency.
        for tt = 1:nSFs
            if (ispref('SpatioSpectralStimulator','SACCMaterials'))
                testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
                testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration',measureDate,dataType,whichChannel);
                testFilename = GetMostRecentFileName(testFiledir,append(num2str(targetCyclePerDeg{tt}),'cpd_crop'));

                % We save all images here. The array looks like {dataType,
                % channel, SF}.
                images{dd,cc,tt} = imread(testFilename);
            else
                error('Cannot find data file');
            end
        end
    end
end

%% Plot the camera images.
PLOTIMAGE = false;
if (PLOTIMAGE)
    for dd = 1:nDataTypes
        % We will make two figures.
        figure;
        figurePosition = [0 0 800 800];

        for cc = 1:nChannels
            set(gcf,'position',figurePosition);
            for tt = 1:nSFs
                subplot(5,5,tt + nSFs*(cc-1));
                imshow(images{dd,cc,tt});
                title(sprintf('%d',tt+nSFs*(cc-1)))
            end
        end
    end
end

%% Plot the sliced images.
PLOTSLICEDIMAGE = true;
if (PLOTSLICEDIMAGE)
    for dd = 1:nDataTypes
        for cc = 1:nChannels
            % Make a new figure per each channel.
            figure;
            figurePosition = [0 0 800 800];
            set(gcf,'position',figurePosition);

            % Add a grand title of the figure.
            sgtitle(sprintf('%s - %s (%s nm)',dataTypeSetting{dd},channels{cc},peakWls{cc}),'FontSize',15);

            % We will set the min peak distance differently to pick the peaks
            % correct for contrast calculation.
            if (dd == 1 && cc == 3)
                minPeakDistance = [20, 22, 10, 5, 4];
            elseif(dd == 2 && cc == 4)
                minPeakDistance = [25, 15, 10, 5, 4];
            else
                minPeakDistance = [20, 15, 10, 5, 4];
            end

            % Calculate the contrasts here.
            for tt = 1:nSFs
                subplot(5,1,tt);
                contrastsRawTemp = GetImgContrast(images{dd,cc,tt},'minPeakDistance',minPeakDistance(tt));
                contrastsRaw{dd,cc,tt} = contrastsRawTemp;
                meanContrasts{dd,cc,tt} = mean(contrastsRawTemp);
                stdErrorContrasts{dd,cc,tt} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
            end

            % Print out progress.
            fprintf('Progress - (%d/%d) \n',cc+nChannels*(dd-1),nChannels*nDataTypes);

            % Save the plot if you want.
            SAVETHEPLOT = false;
            if (SAVETHEPLOT)
                testFileFormat = '.tiff';
                testFilename = sprintf('%s_%s',dataTypeSetting{dd},channels{cc});
                saveas(gcf,append(testFilename,testFileFormat));
                disp('Plot has been saved successfully!');
            end
        end
    end
end

%% Plot the contrasts results - Raw.
%
% Make a new figure.
figure; hold on;
SFs = cell2mat(targetCyclePerDeg);
xticks(SFs);
xticklabels(SFs);
xlabel('Spatial Frequency (cpd)');
ylabel('Mean Contrast');
title('Without optical system (Raw)')
ylim([0.5 1.05]);

% Raw image
colorLines = {'b','c','g',[0.8 0.6 0],'r'};
for cc = 1:nChannels
    numDataType = 1;
    meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
    stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
    plot(SFs,meanContrastTemp,'color',colorLines{cc});
    errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
end
r = get(gca,'Children');
set(r([1:2:7]),'LineStyle','none');
legend(flip(r([2:2:10])),append(peakWls,' nm'),'location','southwest');

%% Plot the contrasts results - SACCSFA.
%
% Make a new figure.
figure; hold on;
SFs = cell2mat(targetCyclePerDeg);
xticks(SFs);
xticklabels(SFs);
xlabel('Spatial Frequency (cpd)');
ylabel('Mean Contrast');
title('With optical system (SACCSFA)')
ylim([0.5 1.05]);

for cc = 1:nChannels
    numDataType = 2;
    meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
    stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
    plot(SFs,meanContrastTemp,'color',colorLines{cc});
    errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
end

s = get(gca,'Children');
set(s([1:2:7]),'LineStyle','none');
legend(flip(s([2:2:10])),append(peakWls,' nm'),'location','southwest');

%% Plot the contrasts comparing within the same channel.
%
% Make a new figure.
figure; hold on;
SFs = cell2mat(targetCyclePerDeg);
lineStyles = {'-','--'};

for cc = 1:nChannels
    subplot(2,3,cc); hold on;

    for dd = 1:nDataTypes
        numDataType = dd;
        meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
        stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
        plot(SFs,meanContrastTemp,lineStyles{dd},'color',colorLines{cc},'linewidth',1.2);
        errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
    end
    % Add plot details.
    xticks(SFs);
    xticklabels(SFs);
    xlabel('Spatial Frequency (cpd)');
    ylabel('Mean Contrast');
    ylim([0.5 1.05]);
    title(append(peakWls(cc),' nm'));

    % Legend.
    ss = get(gca,'Children');
    set(ss([1,3]),'LineStyle','none');
    legend(flip(ss([2,4])),dataTypeSetting,'location','southwest');
end

%% Plot spectrum used.
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
    testFiledir = fullfile(testFiledir,'Calibration');
    testFilename = GetMostRecentFileName(testFiledir,'SACCPrimary1');

    % We save all images here. The array looks like {dataType,
    % channel, SF}.
    data = load(testFilename);
end

calData = data.cals{end};
S = calData.rawData.S;
wls = SToWls(S);
spds = calData.processedData.P_device;
spdsUsed = spds(:,[2 5 9 12 15]);

% Plot it.
figure; hold on;
plot(wls,spds,'k-','linewidth',0.8);
plot(wls,spdsUsed,'-','color',[1 0 0 0.3],'linewidth',4);
xlabel('Wavelength (nm)');
ylabel('Spectral power');
xticks([380:80:780]);
xticklabels([380:80:780]);
ylim([0 max(max(spds))*1.05])
f = get(gca, 'children');
legend(f(flip([1 17])),'All channels','Used channels')
title('Channels used for measuring chromatic aberration');
