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
%    11/22/23   smo    - Cleared up a lot and now it is working by
%                        calculating the MTF of both camera and SACCSFA
%                        within this routine.

%% Initialize.
clear; close all;

%% Set variables.
%
% Set spatial frequency levels.
targetCyclePerDeg = {3,6,9,12,18};
nSFs = length(targetCyclePerDeg);

% Set the index for searching for the channel to compare MTF with the SACCSFA.
% The peak wavelengths of the SACCSFA were 422, 476, 530, 592, 658 nm at .
peaks_spd_SACCSFA = [422 476 530 592 658];
idxChannelTarget = [2 3 5 6 8];
nChannelsTarget = length(idxChannelTarget);

% Choose which contrast calculation method to use. 
%
% Set 'optionContrastCalMethod' to 1 will show the results using the
% contrasts calculated directly from the intensity profile of each image.
%
% Set it 2 shows the results by doing sine fitting to the intensity profile
% and calculating the contrasts from the fitting.
optionContrastCalMethod = 1;
switch optionContrastCalMethod
    case 1
        contrastCalMethod = 'MeanIntensityProfile';
    case 2
        contrastCalMethod = 'Sinefit';
end

% Set additional analysis options.
DoFourierTransform = false;
plotIntensityProfile = false;

%% Get the peak wavelength of the Combi-LED.
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
peaks_spd_camera = FindPeakSpds(spd,'verbose',false);

%% 1) Calculate the MTF (SACCSFA).
%
% Set viewing media to load the images.
viewingMedia = 'SACCSFA';

% Load all images here.
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

% Get the available channels by getting the folder names.
countChannel = 1;
for cc = 1:length(channelFolderList)
    channelFoldernameTemp = channelFolderList(cc).name;
    
    % Extract the number of channels only.
    folderNamePattern = 'Ch';
    if strncmp(channelFoldernameTemp,folderNamePattern,length(folderNamePattern))
        numChannels(countChannel) = str2num(cell2mat(regexp(channelFoldernameTemp, '\d+', 'match')));
        channelOptions{countChannel} = channelFoldernameTemp;
        countChannel = countChannel+1;
    end
end

% Sort the channel options in an ascending order.
[numChannelsSorted I] = sort(numChannels,'ascend');

% Sort the channel options in a ascending order here.
channelOptions = channelOptions(I);

% Load all images here for all channels and spatial frequencies.
for cc = 1:nChannelsTarget
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % Get the images of all spatial frequency.
    for ss = 1:nSFs
        testFilenameTemp = GetMostRecentFileName(oneChannelFileDir,...
            append(num2str(targetCyclePerDeg{ss}),'cpd_crop'));
        
        % We save all images here.
        images{ss} = imread(testFilenameTemp);
        
        % Set min distance between adjacent peaks.
        if ss == 1
            minPeakDistance = 35;
        elseif ss == 2
            minPeakDistance = 20;
        else
            minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        if (plotIntensityProfile)
            subplot(nSFs,1,ss);
            title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
        end
        
        % Calculate contrasts.
        [contrastsTemp IP{cc,ss}] = GetImgContrast(images{ss},'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
        meanContrastsOneChannel(ss) = mean(contrastsTemp);
    end
    
    % Collect the mean contrast results.
    meanContrasts_SACCSFA(cc,:) = meanContrastsOneChannel;
end

%% Fit sine function to the signal (SACCSFA).
if strcmp(viewingMedia,'Print')
    IP = IP(idxChannelTarget,:);
end

% Fit sine signal here.
for ss = 1:nSFs
    for cc = 1:nChannelsTarget
        % Set initial frequency for fitting sine wave.
        switch ss
            case 1
                f0Options = [2.63158 2.10526 2.10526 2.42105 2.52632];
            case 2
                f0Options = [7.68421 7.26316 7.05263 7.26316 6.63158];
            case 3
                f0Options = [12.73684 12.52632 12.94747 11.5 13.89474];
            case 4
                f0Options = [16 16.94737 19.15789 18.36842 17.73684];
            case 5
                f0Options = [29.52632 30.36842 30.42105 29.31579 29];
        end
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Fit happens here.
        signalToFit = IP{cc,ss};
        [params_SACSSFA{cc,ss}, fittedSignal_SACCSFA{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
end

% Plot the results.
for ss = 1:nSFs
    % Make a new figure per each spatial frequency.
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('%d cpd (%s)',targetCyclePerDeg{ss},viewingMedia));
    
    % Loop over the channels.
    for cc = 1:nChannelsTarget
        subplot(5,1,cc); hold on;
        title(sprintf('%d nm',peaks_spd_SACCSFA((cc))));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(IP{cc,ss},'b-');
        
        % Fitted signal.
        plot(fittedSignal_SACCSFA{cc,ss},'r-');
        legend('Origianl','Fit');
    end
end

% Calculate contrast from the sine fitted curve.
for ss = 1:nSFs
    for cc = 1:nChannelsTarget
        paramsTemp = params_SACSSFA{cc,ss};
        A = paramsTemp(1);
        B = paramsTemp(4);
        contrast = A/B;
        if contrast > 1
            contrast = 1;
        elseif contrast < 0
            contrast = 0;
        end
        contrastsFit_SACCSFA(cc,ss) = contrast;
    end
end

%% 2) Calculate the MTF (Camera).
%
% Set the viewing media for the camera MTF measurement. We used the printed
% target so the images were saved in the folder 'Print'.
viewingMedia = 'Print';

% Load all images.
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

% Get the available channels by getting the folder names.
countChannel = 1;
for cc = 1:length(channelFolderList)
    channelFoldernameTemp = channelFolderList(cc).name;
    
    % Extract the number of channels only.
    folderNamePattern = 'Ch';
    if strncmp(channelFoldernameTemp,folderNamePattern,length(folderNamePattern))
        numChannels(countChannel) = str2num(cell2mat(regexp(channelFoldernameTemp, '\d+', 'match')));
        channelOptions{countChannel} = channelFoldernameTemp;
        countChannel = countChannel+1;
    end
end

% Sort the channel options in an ascending order.
[numChannelsSorted I] = sort(numChannels,'ascend');

% Sort the channel options in a ascending order here.
channelOptions = channelOptions(I);

% Load all images here for all channels and spatial frequencies.
for cc = 1:nChannels
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % Get the images of all spatial frequency.
    for ss = 1:nSFs
        testFilenameTemp = GetMostRecentFileName(oneChannelFileDir,...
            append(num2str(targetCyclePerDeg{ss}),'cpd_crop'));
        
        % We save all images here.
        images{ss} = imread(testFilenameTemp);
        
        % Set min distance between adjacent peaks.
        if ss == 1
            minPeakDistance = 35;
        elseif ss == 2
            minPeakDistance = 20;
        else
            minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        if (plotIntensityProfile)
            subplot(nSFs,1,ss);
            title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
        end
        
        % Calculate contrasts.
        [contrastsTemp IP{cc,ss}] = GetImgContrast(images{ss},'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
        meanContrastsOneChannel(ss) = mean(contrastsTemp);
    end
    
    % Collect the mean contrast results.
    meanContrasts_camera(cc,:) = meanContrastsOneChannel;
end

%% Fit sine function to the signal (Camera).
if strcmp(viewingMedia,'Print')
    IP = IP(idxChannelTarget,:);
end

% Fit sine signal here.
%
% Loop over the spatial frequency.
for ss = 1:nSFs
    % Loop over the channels.
    for cc = 1:nChannelsTarget
        
        % Set initial frequency for fitting sine wave. Fitting results are
        % extremely sensitive how we set the initial guess of its frequency. We
        % recommend setting it close to its fundamental frequency.
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
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Fit happens here.
        signalToFit = IP{cc,ss};
        [params_camera{cc,ss}, fittedSignal_camera{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
    
    % Show the progress.
    fprintf('Fitting progress - (%d/%d) \n',ss,nSFs);
end

% Plot the results.
for ss = 1:nSFs
    % Make a new figure per each spatial frequency.
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('%d cpd (%s)',targetCyclePerDeg{ss},viewingMedia));
    
    % Loop over the channels.
    for cc = 1:nChannelsTarget
        subplot(5,1,cc); hold on;
        title(sprintf('%d nm',peaks_spd_camera(idxChannelTarget(cc))));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(IP{cc,ss},'b-');
        
        % Fitted signal.
        plot(fittedSignal_camera{cc,ss},'r-');
        legend('Origianl','Fit');
    end
end

% Calculate contrast from the sine fitted curve.
for ss = 1:nSFs
    for cc = 1:nChannelsTarget
        paramsTemp = params_camera{cc,ss};
        A = paramsTemp(1);
        B = paramsTemp(4);
        contrast = A/B;
        if contrast > 1
            contrast = 1;
        elseif contrast < 0
            contrast = 0;
        end
        contrastsFit_camera(cc,ss) = contrast;
    end
end

%% Calculate the compensated MTF (Camera).
%
% Here we compensate the limitation of using printed paper by dividing the
% MTF of a single cycle to the camera MTF.
%
% Get the pre-saved MTF of the signle cycle. It contains a total of 8
% elements which match with the number of channels in the Combi-LED
% projector.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration');
testFilename = 'Contrast_SingleCycle_combiLED.mat';
data = load(fullfile(testFiledir,testFilename));
contrast_singleCycle = data.contrastSingleCyclePerChannel';

% Normalize the contrast by dividing the single cycle contrast.
meanContrasts_cameraNorm = meanContrasts_camera./contrast_singleCycle;
contrastsFit_cameraNorm = contrastsFit_camera./contrast_singleCycle(idxChannelTarget);

% Set the contrast within the range.
maxContrast = 1;
meanContrasts_cameraNorm(find(meanContrasts_cameraNorm > maxContrast)) = 1;
contrastsFit_cameraNorm(find(contrastsFit_cameraNorm > maxContrast)) = 1;

%% 3) Calculate the compensated MTF (SACCSFA).
%
% We used two different methods to calculate contrast. Choose either one to
% plot the results. It was chosen at the very beginning of this routine.
switch contrastCalMethod
    case 'MeanIntensityProfile'
        contrast_camera = meanContrasts_cameraNorm(idxChannelTarget,:);
        contrast_SACCSFA = meanContrasts_SACCSFA;
    case 'Sinefit'
        contrast_camera = contrastsFit_cameraNorm;
        contrast_SACCSFA = contrastsFit_SACCSFA;
end

% Make a loop to plot the results of each channel.
figure; clf;
figureSize = [0 0 1000 500];
set(gcf,'position',figureSize);
sgtitle('MTF: Camera vs. SACCSFA', 'fontsize', 15);

for cc = 1:nChannelsTarget
    subplot(2,3,cc); hold on;
    
    % Camera MTF.
    meanContrastsOneChannel = contrast_camera(cc,:);
    plot(cell2mat(targetCyclePerDeg),meanContrastsOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % SACCSFA MTF.
    plot(cell2mat(targetCyclePerDeg),contrast_SACCSFA(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','r','markersize',10);
    
    legend('Raw','SACCSFA','location','southeast','fontsize',11);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(cc)), 'fontsize', 15);
end

% Now plot the MTF of SACCSFA under assumption using a perfect camera.
%
% Here we divide the MTF by camera MTF.
contrast_SACCSFA_compensated = contrast_SACCSFA./contrast_camera;

% Set the contrast within the range.
maxContrast = 1;
contrast_SACCSFA_compensated(find(contrast_SACCSFA_compensated > maxContrast)) = 1;

% Make a new figure.
figure; clf;
set(gcf,'position',figureSize);

% Make a loop to plot the results of each channel.
for cc = 1:length(contrast_SACCSFA_compensated)
    subplot(2,3,cc); hold on;
    meanContrastsOneChannel = contrast_SACCSFA_compensated(cc,:);
    plot(cell2mat(targetCyclePerDeg),meanContrastsOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(cc)), 'fontsize', 15);
    legend('SACCSFA','location','southeast','fontsize',11);
end

%% Transverse Chromatic Aberration (TCA) - (Camera).
%
% Here we will compare intensity profile across the channels.
figure; hold on;
figurePosition = [0 0 1000 1000];
set(gcf,'position',figurePosition);
sgtitle('Check spatial position of the waves over different channels');

% Make a loop to plot all combinations, channel x spatial frequency.
%
% Spatial frequency.
peaks_spd_camera_compare = peaks_spd_camera(idxChannelTarget);
for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannelsTarget
        ESFTemp = IP{cc,ss};
        plot(ESFTemp);
        
        % Generate texts for the legend.
        legendHandlesESF{cc} = append(num2str(peaks_spd_camera_compare(cc)),' nm');
    end
    
    % Set each graph in format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandlesESF,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('Normalized dRGB','fontsize',12);
end

%% For better fitting the sine curve, we can search the initial guess of frequency.
%
% We already found optimal values of initial frequency settings, but we
% keep this routine here for further use. Eventually we may want to put
% this part inside the fitting routine.
%
% Search a value using grid-search.
FINDINITIALFREQUENCYTOFIT = false;

if (FINDINITIALFREQUENCYTOFIT)
    nFits = 20;
    f0_lb = 3;
    f0_ub = 5;
    f0Range = linspace(f0_lb,f0_ub,nFits);
    
    % Set the wave to fit.
    SF = 1;
    originalSignals = IP(:,SF);
    
    for cc = 1:nChannelsTarget
        % Make a new figure per each channel.
        figure; hold on;
        figurePosition = [0 0 1000 1000];
        set(gcf,'position',figurePosition);
        sgtitle(sprintf('%d nm', peaks_spd_camera(idxChannelTarget(cc))));
        
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
