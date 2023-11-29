% SACC_GetImageContrastWithFunction.
%
% It calculates the image contrasts using the function.

% History:
%    06/13/23   smo    - Wrote it.
%    06/27/23   smo    - Now we load all images and analyze it together.
%    07/13/23   smo    - Added a plot to compare the actual peaks of the
%                        spectrums used in each projector.

%% Initialize.
clear; close all;

%% Set variables.
%
% Initial measurements were made on 0613.
targetCyclePerDeg = {3, 6, 9, 12, 18};
projectorSettings = {'SACCSFA'};
measureDate = '0906';
focusedImage = false;

%% Plot the spectrum used.
PLOTSPECTRUM = false;
if (PLOTSPECTRUM)
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Calibration');
        testFilename = GetMostRecentFileName(testFiledir,'SACCPrimary1');
        
        % We save all images here. The array looks like {dataType,
        % channel, SF}.
        data = load(testFilename);
    end
    
    % Here we read out both old and new projector calibration files. For new
    % projector, we read the most recent one. For old projector, we read the
    % last one measured which is stored in 17th of the calibration file in
    % SACCPriamry1.
    idxFileOldProjector = 17;
    calData_oldProjector = data.cals{idxFileOldProjector};
    calData_newProjector = data.cals{end};
    S = calData_newProjector.rawData.S;
    wls = SToWls(S);
    spds_newProjector = calData_newProjector.processedData.P_device;
    spds_oldProjector = calData_oldProjector.processedData.P_device;
    
    % Find the peak wavelength.
    peaks_spd_new = FindPeakSpds(spds_newProjector);
    
    % Get spds of the channels used.
    %
    % We used the same channels for both old and new projector [2 5 9 12 15] on
    % the date of 0613. We will set these numbers differently so that the
    % actual peaks are matched each other.
    switch measureDate
        case '0613'
            numChannelUsed_newProjector = [2 5 9 12 15];
            numChannelUsed_oldProjector = [2 5 9 12 15];
        case '0714'
            numChannelUsed_newProjector = [1 3 7 12 15];
            numChannelUsed_oldProjector = [1 3 7 11 14];
        case '0718'
            numChannelUsed_newProjector = [1 3 7 10 15];
            numChannelUsed_oldProjector = [1 3 7 9 14];
        case '0719'
            numChannelUsed_newProjector = [1 3 7 10 15];
            numChannelUsed_oldProjector = [1 3 7 9 14];
    end
    
    spdsUsed_newProjector = spds_newProjector(:,numChannelUsed_newProjector);
    spdsUsed_oldProjector = spds_oldProjector(:,numChannelUsed_oldProjector);
    
    % Peaks of each spectrums.
    peaksUsed_newProjector = FindPeakSpds(spdsUsed_newProjector,'verbose',false);
    peaksUsed_oldProjector = FindPeakSpds(spdsUsed_oldProjector,'verbose',false);
    
    % Plot it.
    figure;
    figPosition = [0 0 1500 400];
    set(gcf,'position',figPosition);
    sgtitle('Subprimary channels used for measuring chromatic aberration','fontsize',15);
    
    % New projector.
    subplot(1,3,1); hold on;
    plot(wls,spds_newProjector,'k-','linewidth',0.8);
    plot(wls,spdsUsed_newProjector,'-','color',[1 0 0 0.3],'linewidth',4);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    xticklabels([380:80:780]);
    ylim([0 max(max(spds_newProjector))*1.05])
    f = get(gca, 'children');
    legend(f(flip([1 17])),'All channels','Used channels',...
        'location','northwest','fontsize',13)
    title('SACCSFA','fontsize',15);
    subtitle(sprintf('Peaks = (%s) nm',num2str(peaksUsed_newProjector)),'fontsize',14);
    
    % Old projector.
    subplot(1,3,2); hold on;
    plot(wls,spds_oldProjector,'k-','linewidth',0.8);
    plot(wls,spdsUsed_oldProjector,'-','color',[0 1 0 0.3],'linewidth',4);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    xticklabels([380:80:780]);
    ylim([0 max(max(spds_newProjector))*1.05])
    f = get(gca, 'children');
    legend(f(flip([1 17])),'All channels','Used channels',...
        'location','northwest','fontsize',13)
    title('Raw (Old Projector)','fontsize',15);
    subtitle(sprintf('Peaks = (%s) nm',num2str(peaksUsed_oldProjector)),'fontsize',14);
    
    % Comparison New vs. Old projector.
    subplot(1,3,3), hold on;
    plot(wls,spdsUsed_newProjector,'-','color',[1 0 0 0.3],'linewidth',4);
    plot(wls,spdsUsed_oldProjector,'-','color',[0 1 0 0.3],'linewidth',4);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    xticklabels([380:80:780]);
    ylim([0 max(max(spds_newProjector))*1.05])
    f = get(gca, 'children');
    legend(flip(f([1 6])),'SACCSFA','Raw (Old Projector)',...
        'location','northwest','fontsize',13)
    title('SACCSFA vs. Raw','fontsize',15);
    
    % Save the peak wavelengths in string. We will use this for legend in the
    % following plot.
    for pp = 1:length(peaksUsed_newProjector)
        peakWls{pp} = num2str(peaksUsed_newProjector(pp));
    end
else
    peakWls = {'422', '476', '530', '592', '658'};
end

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
        
        for tt = 1:nSFs
            testFiledirTemp = fullfile(testFiledir,channelTemp);
            % You can load separately focused image if you want.
            if (focusedImage)
                try
                    testFilename = GetMostRecentFileName(testFiledirTemp,append(num2str(targetCyclePerDeg{tt}),'cpd_focused_crop'));
                catch
                    % If there is no such file name, just load regualr
                    % image.
                    disp('No such file name found. Regular image file will be loaded');
                    testFilename = GetMostRecentFileName(testFiledirTemp,append(num2str(targetCyclePerDeg{tt}),'cpd_crop'));
                end
            else
                % Get the file name of the images.
                testFilename = GetMostRecentFileName(testFiledirTemp,append(num2str(targetCyclePerDeg{tt}),'cpd_crop'));
            end
            
            % We save all images here. The array looks like {dataType,
            % channel, SF}.
            images{dd,cc,tt} = imread(testFilename);
        end
    end
end

%% Plot the camera images.
PLOTIMAGE = false;
if (PLOTIMAGE)
    for dd = 1:nProjectorSettings
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
    for dd = 1:nProjectorSettings
        for cc = 1:nChannels
            % Make a new figure per each channel.
            figure;
            figurePosition = [0 0 800 800];
            set(gcf,'position',figurePosition);
            
            % Add a grand title of the figure.
            sgtitle(sprintf('%s - (%s nm)',projectorSettings{dd},peakWls{cc}),'FontSize',15);
            
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
                [contrastsRawTemp ESFTemp] = GetImgContrast(images{dd,cc,tt},'minPeakDistance',minPeakDistance(tt));
                contrastsRaw{dd,cc,tt} = contrastsRawTemp;
                meanContrasts{dd,cc,tt} = mean(contrastsRawTemp);
                stdErrorContrasts{dd,cc,tt} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
                ESF{cc,tt} = ESFTemp;
            end
            
            % Print out progress.
            fprintf('Progress - (%d/%d) \n',cc+nChannels*(dd-1),nChannels*nProjectorSettings);
            
            % Save the plot if you want.
            SAVETHEPLOT = false;
            if (SAVETHEPLOT)
                testFileFormat = '.tiff';
                testFilename = sprintf('%s_%s (nm)',projectorSettings{dd},peakWls{cc});
                saveas(gcf,append(testFilename,testFileFormat));
                disp('Plot has been saved successfully!');
            end
        end
    end
end

%% Plot the contrasts results - Raw.
%
% Define color of the lines on the plot.
colorLines = {'b','c','g',[0.8 0.6 0],'r'};

% Make a new figure.
if ismember('Raw',projectorSettings)
    figure; hold on;
    SFs = cell2mat(targetCyclePerDeg);
    xticks(SFs);
    xticklabels(SFs);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrast','fontsize',15);
    title('Raw (Old projector)','fontsize',15)
    ylim([0 1.05]);
    
    % Raw image
    for cc = 1:nChannels
        numDataType = 1;
        meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
        stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
        plot(SFs,meanContrastTemp,'o-','color',colorLines{cc},'linewidth',1.5);
        errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
    end
    r = get(gca,'Children');
    set(r([1:2:7]),'LineStyle','none');
    legend(flip(r([2:2:10])),append(peakWls,' nm'),'location','southwest','fontsize',15);
end

%% Plot the contrasts results - SACCSFA.
%
% Make a new figure.
if ismember('SACCSFA',projectorSettings)
    figure; hold on;
    SFs = cell2mat(targetCyclePerDeg);
    xticks(SFs);
    xticklabels(SFs);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrast','fontsize',15);
    title('SACCSFA','fontsize',15)
    ylim([0 1.05]);
    
    for cc = 1:nChannels
        if length(projectorSettings)==1
            numDataType = 1;
        else
            numDataType = 2;
        end
        meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
        stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
        plot(SFs,meanContrastTemp,'o-','color',colorLines{cc},'linewidth',1.5);
        errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
    end
    
    s = get(gca,'Children');
    set(s([1:2:7]),'LineStyle','none');
    legend(flip(s([2:2:10])),append(peakWls,' nm'),'location','southwest','fontsize',15);
end

%% Plot the contrasts comparing within the same channel.
%
% Make a new figure.
figure; hold on;
SFs = cell2mat(targetCyclePerDeg);
lineStyles = {'o-','o--'};

for cc = 1:nChannels
    subplot(2,3,cc); hold on;
    
    for dd = 1:nProjectorSettings
        numDataType = dd;
        meanContrastTemp = cell2mat(squeeze(meanContrasts(numDataType,cc,:)));
        stdErrorContrastTemp = cell2mat(squeeze(stdErrorContrasts(numDataType,cc,:)));
        plot(SFs,meanContrastTemp,lineStyles{dd},'color',colorLines{cc},'linewidth',1.2);
        errorbar(SFs,meanContrastTemp,stdErrorContrastTemp,'color',colorLines{cc});
    end
    % Add plot details.
    xticks(SFs);
    xticklabels(SFs);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrast','fontsize',15);
    ylim([0 1.05]);
    title(append(peakWls(cc),' nm'),'fontsize',15);
    
    % Add legend.
    ss = get(gca,'Children');
    if length(projectorSettings) == 1
        set(ss(1),'LineStyle','none');
    else
        set(ss([1,3]),'LineStyle','none');
        legend(flip(ss([2,4])),projectorSettings,'location','southwest','fontsize',15);
    end
    
    % Collect the mean contrasts data in an array.
    meanContrasts_all(:,cc) = meanContrastTemp;
end

%% Compare ESF across the channels.
%
% The size of array 'ESF' is 8 (channels) x 5 (spatial frequency).
figure; hold on;
figurePosition = [0 0 1000 1000];
set(gcf,'position',figurePosition);
sgtitle('Check spatial position of the waves over different channels');

% Make a loop to plot all combinations, channel x spatial frequency.
%
% Spatial frequency.
for dd = 1:nSFs
    subplot(5,1,dd); hold on;
    
    % Channel.
    for cc = 1:nChannels
        ESFTemp = ESF{cc,dd};
        plot(ESFTemp);
        
        % Generate texts for the legend.
        legendHandlesESF{cc} = append(num2str(peakWls{cc}),' nm');
    end
    
    % Set each graph in format.
    title(sprintf('%d cpd',targetCyclePerDeg{dd}),'fontsize',15);
    legend(legendHandlesESF,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('Normalized dRGB','fontsize',12);
end

%% Fit sine function to the signal.
%
% Show FFT results if you want.
DoFourierTransform = false;

% Fit sine signal here.
%
% Loop over the spatial frequency.
for dd = 1:nSFs
    % Loop over the channels.
    for cc = 1:nChannels
        
        % Set initial frequency for fitting sine wave. Fitting results are
        % extremely sensitive how we set the initial guess of its frequency. We
        % recommend setting it close to its fundamental frequency.
        switch dd
            case 1
                f0Options = [2.63158 2.10526 2.10526 2.42105 2.52632];
            case 2
                f0Options = [7.68421 7.26316 7.05263 7.26316 6.63158];
            case 3
                f0Options = [12.73684 12.52632 12.94747 11.5 13.89474];
            case 4
                f0Options = [16 16.94737   19.15789  18.36842  17.73684];
            case 5
                f0Options = [29.52632 30.36842 30.42105 29.31579 29];
        end
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Get one signal to fit.
        signalTemp = ESF{cc,dd};
        
        % Fit happens here.
        [params{cc,dd}, fittedSignal{cc,dd}] = FitSineWave(signalTemp,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
    fprintf('Fitting progress - (%d/%d) \n',dd,nSFs);
end

% Plot the results.
for dd = 1:nSFs
    % Make a new figure per each spatial frequency.
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('%d cpd',targetCyclePerDeg{dd}));
    
    % Loop over the channels.
    for cc = 1:nChannels
        subplot(5,1,cc); hold on;
        title(sprintf('%s nm',peakWls{cc}));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(ESF{cc,dd},'b-');
        
        % Fitted signal.
        plot(fittedSignal{cc,dd},'r-');
        legend('Origianl','Fit');
    end
end

%% we search the initial guess of frequency to fit sine curve.
%
% Search a value using grid-search.
nFits = 20;
f0_lb = 2;
f0_ub = 4;
f0Range = linspace(f0_lb,f0_ub,nFits);

% Set the wave to fit.
SF = 1;
originalSignals = ESF(:,SF);

for cc = 1:nChannels
    % Make a new figure per each channel.
    figure; hold on;
    figurePosition = [0 0 1000 1000];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('%s nm', peakWls{cc}));
    
    for ff = 1:nFits
        f0 = f0Range(ff);
        originalSignal = originalSignals{cc};
        [~, fittedSignalOne] = FitSineWave(originalSignal,'f0',f0,'verbose',false,'FFT',false);
        
        % Original.
        subplot(5,4,ff); hold on;
        plot(originalSignal,'b-');
        plot(fittedSignalOne,'r-');
        title(sprintf('f0 = %.5f',f0));
        legend('Origianl','Fit');
        ylim([0 max(originalSignal)*1.05]);
        
        clear f0;
    end
end
