% SACC_ChromaticAberrationAnalyze.
%
% This is to analyze chromatic aberration of the SACCSFA system. This is
% universal routine that can be used for measured image of SACCSFA, Raw
% projector, Pritned image.
%
% This has been developed based on the code,
% SACC_ContrastPrintedImage_combiLED, which analyzes the image captured
% using the Printed image.
%
% See also:
%    SACC_GetCameraImageContrast, SACC_ContrastPrintedImage_combiLED

% History:
%    11/17/23   smo    - Wrote it to use the routine for all viewing
%                        media (SACCSFA, Print, RawProjector).

%% Initialize.
clear; close all;

%% Set variables.
targetCyclePerDeg = {3,6,9,12,18};
nSFs = length(targetCyclePerDeg);

% Set which viewing media data to analyze.
numViewingMedia = 2;
switch numViewingMedia
    case 1
        viewingMedia = 'SACCSFA';
    case 2
        viewingMedia = 'Print';
    case 3
        viewingMedia = 'RawProjector';
end

%% Get the peak wave length values.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration');
testFilename = 'spd_combiLED.mat';
spdData = load(fullfile(testFiledir,testFilename));

% Extract the spd data and flip left to right, becasue the measurement was
% done from Ch8 (high, 652 nm) to Ch1 (low, 406 nm).
spd = spdData.spd;
spd = fliplr(spd);
nChannels = size(spd,2);
S = [380 2 201];
wls = SToWls(S);
peaks_spd = FindPeakSpds(spd,'verbose',false);

%% Load SACCSFA MTF results to compare.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration');
testFilename = 'MTF_SACCSFA.mat';
data_SACCSFA = load(fullfile(testFiledir,testFilename));
meanContrasts_SACCSFA = data_SACCSFA.meanContrasts_all;

%% Load all images here.
%
% Get availble filelist.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration',viewingMedia);
folderList = dir(testFiledir);
folderList = folderList(4:end);

% Find the most recent measurement.
for i = 1:numel(folderList)
    folderDate(i) = datetime(folderList(i).name);
end
[recentDate, idxRecentDate] = max(folderDate);
recentFolderName = folderList(idxRecentDate).name;
recentTestFiledir = fullfile(testFiledir,recentFolderName);

% Find available channels.
channelFolderList = dir(recentTestFiledir);
channelFolderList = channelFolderList(4:end);

nChannels = numel(channelFolderList);
nSFs = length(targetCyclePerDeg);

% Load all images here for all channels and spatial frequencies.
for cc = 1:nChannels
    % Set the folder of each channel.
    channelOptions{cc} = channelFolderList(cc).name;
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % Loop over for all spatial frequency.
    for ss = 1:nSFs
        testFilenameTemp = GetMostRecentFileName(oneChannelFileDir,...
            append(num2str(targetCyclePerDeg{ss}),'cpd_crop'));
        
        % We save all images here.
        images{ss} = imread(testFilenameTemp);
    end
    
    %% Plot the sliced images.
    %
    % Contrast calculation happens here.
    for ss = 1:nSFs
        % Set min distance between adjacent peaks.
        if ss == 1
            minPeakDistance = 35;
        elseif ss == 2
            minPeakDistance = 20;
        else
            minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        subplot(nSFs,1,ss);
        title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
        
        % Calculate contrasts.
        [contrastsRawTemp ESF{cc,ss}] = GetImgContrast(images{ss},'minPeakDistance',minPeakDistance);
        contrastsRaw{ss} = contrastsRawTemp;
        meanContrasts{ss} = mean(contrastsRawTemp);
        stdErrorContrasts{ss} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
    end
    
    % Collect the mean contrast results.
    meanContrasts_all(:,cc) = cell2mat(meanContrasts);
end

%% Plot MTF comparing with SACCSFA results.
%
% This is one cycle contrast results. We will calculate the compensated
% contrasts by dividing when plotting the results.
%
% These are the contrats when measuring only one cycle (so, half black on
% the left and the other half as white on the right). Measured on
% 09/05/2023.
contrastSingleCyclePerChannel = [0.8811 0.8756 0.8938 0.8891 0.8863 0.8870 0.8853 0.8868];
if strcmp(viewingMedia,'SACCSFA')
    idxChComparison = [2 3 5 6 8];
    contrastSingleCyclePerChannel = contrastSingleCyclePerChannel(idxChComparison);
end

figure; clf;
figureSize = [0 0 1000 500];
set(gcf,'position',figureSize);
sgtitle('MTF comparison: Camera vs. SACCSFA', 'fontsize', 15);

% Normalize the contrast by dividing the single cycle contrast.
meanContrastsNorm_all = meanContrasts_all./contrastSingleCyclePerChannel;

% Set the contrast within the range.
maxContrast = 1;
meanContrastsNorm_all(find(meanContrastsNorm_all > maxContrast)) = 1;

% Set the index for searching for the channel to compare MTF with the SACCSFA.
% The peak wavelengths of the SACCSFA were 422, 476, 530, 592, 658 nm at .
peaks_spd_SACCSFA = [422 476 530 592 658];
idxChComparison = [2 3 5 6 8];
nChComparison = length(idxChComparison);

% Make a loop to plot the results of each channel.
for cc = 1:nChComparison
    subplot(2,3,cc); hold on;
    
    % Printed pattern.
    meanContrastsTemp = meanContrastsNorm_all(:,cc);
    plot(cell2mat(targetCyclePerDeg),meanContrastsTemp,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % SACCSFA.
    plot(cell2mat(targetCyclePerDeg),meanContrasts_SACCSFA(:,ss),...
        'ko-','markeredgecolor','k','markerfacecolor','r','markersize',10);
    
    legend('Raw','SACCSFA','location','southeast','fontsize',11);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(cc)), 'fontsize', 15);
end

%% Plot the MTF of SACCSFA under assumption using a perfect camera.
%
% Here we divide the MTF by camera MTF.
meanContrastsSACCSFAPerfect = meanContrasts_SACCSFA./meanContrastsNorm_all(idxChComparison);

% Set the contrast within the range.
maxContrast = 1;
meanContrastsSACCSFAPerfect(find(meanContrastsSACCSFAPerfect > maxContrast)) = 1;

% Make a new figure.
figure; clf;
set(gcf,'position',figureSize);

% Make a loop to plot the results of each channel.
for cc = 1:length(meanContrastsSACCSFAPerfect)
    subplot(2,3,cc); hold on;
    meanContrastsTemp = meanContrastsSACCSFAPerfect(:,cc);
    plot(cell2mat(targetCyclePerDeg),meanContrastsTemp,...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(cc)), 'fontsize', 15);
    legend('SACCSFA','location','southeast','fontsize',11);
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
peaks_spd_comparison = peaks_spd(idxChComparison);
for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChComparison
        ESFTemp = ESF{cc,ss};
        plot(ESFTemp);
        
        % Generate texts for the legend.
        legendHandlesESF{cc} = append(num2str(peaks_spd_comparison(cc)),' nm');
    end
    
    % Set each graph in format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
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
for ss = 1:nSFs
    % Loop over the channels.
    for cc = 1:nChComparison
        
        % Set initial frequency for fitting sine wave. Fitting results are
        % extremely sensitive how we set the initial guess of its frequency. We
        % recommend setting it close to its fundamental frequency.
        if strcmp(viewingMedia,'Print')
            switch ss
                case 1
                    f0Options = [3.3684 3.6316 3.6316 3.6316 3.6316];
                case 2
                    f0Options = [7.3333 6 7.1282 7.3333 7.5789];
                case 3
                    f0Options = [13.5641 14.4211 14.8276 14.4211 14.4211];
                case 4
                    f0Options = [29 18.6316 18.6316 18.9474 18];
                case 5
                    f0Options = [29.5 30.01 30.5 29.9474 30.02];
            end
        elseif strcmp(viewingMedia,'SACCSFA')
            switch ss
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
        end
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Get one signal to fit.
        signalTemp = ESF{cc,ss};
        
        % Fit happens here.
        [params{cc,ss}, fittedSignal{cc,ss}] = FitSineWave(signalTemp,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
    fprintf('Fitting progress - (%d/%d) \n',ss,nSFs);
end

% Plot the results.
for ss = 1:nSFs
    % Make a new figure per each spatial frequency.
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('%d cpd',targetCyclePerDeg{ss}));
    
    % Loop over the channels.
    for cc = 1:nChComparison
        subplot(5,1,cc); hold on;
        title(sprintf('%d nm',peaks_spd(idxChComparison(cc))));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(ESF{cc,ss},'b-');
        
        % Fitted signal.
        plot(fittedSignal{cc,ss},'r-');
        legend('Origianl','Fit');
    end
end

%% Here we search the initial guess of frequency to fit sine curve.
%
% Search a value using grid-search.
FINDINITIALFREQUENCYTOFIT = true;

if (FINDINITIALFREQUENCYTOFIT)
    nFits = 20;
    f0_lb = 3;
    f0_ub = 5;
    f0Range = linspace(f0_lb,f0_ub,nFits);
    
    % Set the wave to fit.
    SF = 1;
    originalSignals = ESF(:,SF);
    
    for cc = 1:nChComparison
        % Make a new figure per each channel.
        figure; hold on;
        figurePosition = [0 0 1000 1000];
        set(gcf,'position',figurePosition);
        sgtitle(sprintf('%d nm', peaks_spd(idxChComparison(cc))));
        
        for ff = 1:nFits
            f0 = f0Range(ff);
            originalSignal = originalSignals{cc};
            [~, fittedSignalOne] = FitSineWave(originalSignal,'f0',f0,'verbose',false,'FFT',false);
            
            % Original.
            subplot(5,4,ff); hold on;
            plot(originalSignal,'b-');
            plot(fittedSignalOne,'r-');
            title(sprintf('f0 = %.4f',f0));
            legend('Origianl','Fit');
            ylim([0 max(originalSignal)*1.05]);
            
            clear f0;
        end
    end
end
