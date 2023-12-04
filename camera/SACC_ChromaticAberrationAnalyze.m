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

%% Get the peak wavelengths (SACCSFA).
%
% Load the calibration data. We will load the most recent calibration
% results.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Calibration');
testFilename = 'SACCPrimary1.mat';
calData = load(fullfile(testFiledir,testFilename));
recentCalData = calData.cals{end};

% Get the peaks from the spds.
spd_SACCSFA = recentCalData.processedData.P_device;
peaks_spd_SACCSFA = FindPeakSpds(spd_SACCSFA,'verbose',false);

%% Get the peak wavelength of the Combi-LED.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration');
testFilename = 'spd_combiLED.mat';
spdData = load(fullfile(testFiledir,testFilename));

% Extract the spd data and flip left to right, becasue the measurement was
% done from Ch8 (high, 652 nm) to Ch1 (low, 406 nm).
spd_camera = spdData.spd;
spd_camera = fliplr(spd_camera);
nChannels = size(spd_camera,2);

S = [380 2 201];
wls = SToWls(S);
peaks_spd_camera = FindPeakSpds(spd_camera,'verbose',false);

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
peaks_spd_SACCSFA = sort(peaks_spd_SACCSFA,'ascend');

% Sort the channel options in a ascending order here.
channelOptions = channelOptions(I);

% Load all images here for all channels and spatial frequencies.
nChannels = length(channelOptions);
for cc = 1:nChannels
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % We collect the channel index here.
    idxChannels_SACCSFA(cc) = str2num(cell2mat(regexp(channelOptions{cc},'\d+','match')));
    
    % Make a new figure if we plot the intensity profile.
    if (plotIntensityProfile)
        figure;
        sgtitle(sprintf('%d nm (%s)',peaks_spd_SACCSFA(idxChannels_SACCSFA(cc)),viewingMedia),'fontsize',15);
    end

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
        [contrastsTemp, IP_SACCSFA{cc,ss}] = GetImgContrast(images{ss},'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
        meanContrastsOneChannel(ss) = mean(contrastsTemp);
    end
    
    % Collect the mean contrast results.
    meanContrasts_SACCSFA(cc,:) = meanContrastsOneChannel;
end

%% Fit sine function to the signal (SACCSFA).
%
% Fit sine signal here.
for ss = 1:nSFs
    for cc = 1:nChannels
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
                f0Options = [30.9211 30.36842 30.42105 29.31579 29];
        end
        
        % This is temp options.
        f0Options = ones(1, nChannels);
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Fit happens here.
        signalToFit = IP_SACCSFA{cc,ss};
        [params_SACCSFA{cc,ss}, fittedSignal_SACCSFA{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
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
    for cc = 1:nChannels
        subplot(round(nChannels/2),2,cc); hold on;
        title(sprintf('%d nm',peaks_spd_SACCSFA((cc))));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(IP_SACCSFA{cc,ss},'b-');
        
        % Fitted signal.
        plot(fittedSignal_SACCSFA{cc,ss},'r-');
        legend('Origianl','Fit');
    end
end

% Calculate contrast from the sine fitted curve.
for ss = 1:nSFs
    for cc = 1:nChannels
        paramsTemp = params_SACCSFA{cc,ss};
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
if exist('numChannels')
   clear numChannels;
end
if exist('channelOptions')
   clear channelOptions;
end
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
peaks_spd_camera = sort(peaks_spd_camera,'ascend');

% Sort the channel options in a ascending order here.
channelOptions = channelOptions(I);

% Load all images here for all channels and spatial frequencies.
nChannels = length(channelOptions);
for cc = 1:nChannels
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % Make a new figure if we plot the intensity profile.
    if (plotIntensityProfile)
        figure;
        sgtitle(sprintf('%d nm (%s)',peaks_spd_camera(cc),viewingMedia),'fontsize',15);
    end
    
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
        [contrastsTemp IP_camera{cc,ss}] = GetImgContrast(images{ss},'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
        meanContrastsOneChannel(ss) = mean(contrastsTemp);
    end
    
    % Collect the mean contrast results.
    meanContrasts_camera(cc,:) = meanContrastsOneChannel;
    
    %% Here we get 1 cpd image and calculate contrast to compensate the camera MTF.
    testFilename = '1cpd_crop';
    testFilename = GetMostRecentFileName(oneChannelFileDir,testFilename);
    
    % We save all images here.
    image = imread(testFilename);
    
    % Set min distance between adjacent peaks.
    minPeakDistance = 30;
    
    % Calculate contrasts.
    if (plotIntensityProfile)
        figure;
        title(sprintf('%d nm (%s)',peaks_spd_camera(cc),viewingMedia),'fontsize',15);
        subtitle('1 cpd','fontsize',13);
    end
    [contrastsTemp, IP_camera_1cpd{cc}] = GetImgContrast(image,'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
    meanContrasts_camera_1cpd(cc) = mean(contrastsTemp);
    
    % Now fit sine curve to the 1 cpd to calculate contrast.
    f0 = 10;
    signalToFit = IP_camera_1cpd{cc};
    [params_camera_1cpd{cc}, fittedSignal_camera_1cpd{cc}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
    
    % Clear the initial guess of frequency for next fit.
    clear f0;
    
    % Calculate contrast.
    paramsTemp = params_camera_1cpd{cc};
    A = paramsTemp(1);
    B = paramsTemp(4);
    contrast = A/B;
    contrastsFit_camera_1cpd(cc) = contrast;
end

%% Fit sine function to the signal (Camera).
%
% Fit sine signal here.
%
% Loop over the spatial frequency.
for ss = 1:nSFs
    % Loop over the channels.
    for cc = 1:nChannels
        
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
        
        % This is temp options.
        f0Options = ones(1, nChannels);
        
        % Update initial guess of frequency here.
        f0 = f0Options(cc);
        
        % Fit happens here.
        signalToFit = IP_camera{cc,ss};
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
    for cc = 1:nChannels
        subplot(round(nChannels/2),2,cc); hold on;
        title(sprintf('%d nm',peaks_spd_camera(cc)));
        xlabel('Pixel position');
        ylabel('dRGB');
        ylim([0 220]);
        
        % Original.
        plot(IP_camera{cc,ss},'b-');
        
        % Fitted signal.
        plot(fittedSignal_camera{cc,ss},'r-');
        legend('Origianl','Fit');
    end
end

% Calculate contrast from the sine fitted curve.
for ss = 1:nSFs
    for cc = 1:nChannels
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

%% 2) Plot the raw MTF and compensate it (Camera).
%
% Choose which way to calculate the contrast.
switch contrastCalMethod
    case 'MeanIntensityProfile'
        contrastRaw_camera = meanContrasts_camera;
        contrast_camera_1cpd = meanContrasts_camera_1cpd;
    case 'Sinefit'
        contrastRaw_camera = contrastsFit_camera;
        contrast_camera_1cpd = contrastsFit_camera_1cpd;
end

% Plot the raw camera MTF results.
figure; clf;
figureSize = [0 0 1000 500];
set(gcf,'position',figureSize);
sgtitle('Raw camera MTF', 'fontsize', 15);

for cc = 1:nChannels
    subplot(2,4,cc); hold on;
    
    % Camera MTF.
    contrastRawOneChannel = contrastRaw_camera(cc,:);
    plot(cell2mat(targetCyclePerDeg),contrastRawOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    legend('Raw','location','southeast','fontsize',11);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_camera(cc)), 'fontsize', 15);
end

% Calculate the compensated MTF (Camera).
%
% Here we compensate the limitation of using printed paper by dividing the
% MTF of a 1 cpd to the raw camera MTF. We will calculate the contrasts
% using two different methods.
%
% Normalize the contrast by dividing the single cycle contrast.
meanContrasts_cameraNorm = meanContrasts_camera./contrast_camera_1cpd';
contrastsFit_cameraNorm = contrastsFit_camera./contrastsFit_camera_1cpd';

%% 3) Calculate the compensated MTF (SACCSFA).
%
% We used two different methods to calculate contrast. Choose either one to
% plot the results. It was chosen at the very beginning of this routine.
switch contrastCalMethod
    case 'MeanIntensityProfile'
        contrast_camera = meanContrasts_cameraNorm;
        contrast_SACCSFA = meanContrasts_SACCSFA;
    case 'Sinefit'
        contrast_camera = contrastsFit_cameraNorm;
        contrast_SACCSFA = contrastsFit_SACCSFA;
end

% Here we choose which channel of the combi-LED to compare to each channel
% of SACCSFA. We will choose the one with the closest peak wavelength
% between combi-LED and SACCSFA to compare.
peaks_spd_SACCSFA_test = peaks_spd_SACCSFA(idxChannels_SACCSFA);
nChannelsTest = length(peaks_spd_SACCSFA_test);
for tt = 1:nChannelsTest
    % Get one channel from SACCSFA to find the corresponding channel within
    % the combi-LED.
    peak_spd_SACCSFA = peaks_spd_SACCSFA_test(tt);
    
    % Here we find the channel within the combi-LED that has the minimum
    % difference of the above channel in the SACCSFA.
    [~, idx] = min(abs(peaks_spd_camera - peak_spd_SACCSFA));
    
    % Save out the peak wavelength within the combi-LED that has the
    % closest value to the SACCSFA channel wavelength.
    idx_camera_test(tt) = idx;
    peaks_spd_camera_test(tt) = peaks_spd_camera(idx);
end 

% Make a loop to plot the results of each channel.
figure; clf;
figureSize = [0 0 1000 500];
set(gcf,'position',figureSize);
sgtitle('MTF: Camera vs. SACCSFA', 'fontsize', 15);

for cc = 1:nChannelsTest
    subplot(2,round(nChannelsTest)/2,cc); hold on;
    
    % SACCSFA MTF.
    plot(cell2mat(targetCyclePerDeg),contrast_SACCSFA(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','r','markersize',10);
    
    % Camera MTF.
    %
    % We saved the index of corresponding channel to compare within the
    % combi-LED in 'idx_camera_test'.
    meanContrastsOneChannel = contrast_camera(idx_camera_test(cc),:);
    plot(cell2mat(targetCyclePerDeg),meanContrastsOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',12);
    
    % Save out camera MTF to test.
    %
    % We will use this to compensate the SACCSFA MTF later on.
    contrast_camera_test(cc,:) = meanContrastsOneChannel;
    
    % Plot stuffs.
    legend(sprintf('SACCSFA (%d nm)',peaks_spd_SACCSFA(idxChannels_SACCSFA(cc))),...
        sprintf('Camera (%d nm)',peaks_spd_camera_test(cc)),'location','southeast','fontsize',8);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(idxChannels_SACCSFA(cc))), 'fontsize', 15);
end

% Now plot the MTF of SACCSFA under assumption using a perfect camera.
%
% Here we divide the MTF by camera MTF.
contrast_SACCSFA_compensated = contrast_SACCSFA./contrast_camera_test;

% Set the contrast within the range.
maxContrast = 1;
contrast_SACCSFA_compensated(find(contrast_SACCSFA_compensated > maxContrast)) = 1;

% Make a new figure.
figure; clf;
set(gcf,'position',figureSize);

% Make a loop to plot the results of each channel.
for cc = 1:nChannelsTest
    subplot(2,round(nChannelsTest)/2,cc); hold on;
    meanContrastsOneChannel = contrast_SACCSFA_compensated(cc,:);
    plot(cell2mat(targetCyclePerDeg),meanContrastsOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(idxChannels_SACCSFA(cc))), 'fontsize', 15);
    legend('SACCSFA','location','southeast','fontsize',11);
end

%% Transverse Chromatic Aberration (TCA) - (Camera).
%
% 1) Plot raw intensity profiles.
figure; hold on;
figurePosition = [0 0 1000 1000];
set(gcf,'position',figurePosition);
sgtitle('Raw intensity profile over the channels (Camera)');
minY = -20;
maxY = 220;

% Get the channels used for comparison.
peaks_spd_camera_compare = peaks_spd_camera(idxChannelTarget);

% Make a loop to plot.
for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannelsTarget
        plot(IP_camera{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_camera_compare(cc)),' nm');
        
        % Extract the fitted parameter, phi, for all channels and spatial
        % frequencies.
        idxParamPhi = 3;
        phi_camera(cc,ss) = params_camera{cc,ss}(idxParamPhi);
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% 2) Plot the sine fitted graphs (Camera).
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Fitted intensity profile over the channels (Camera)');

% Loop over Spatial frequency.
for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Loop over Channel.
    for cc = 1:nChannelsTarget
        plot(fittedSignal_camera{cc,ss});
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% 3) Plot the comparison of the parameter phi over the channels.
%
% Define the x-ticks for the plot.
xticksPlot = linspace(1,nChannelsTarget,nChannelsTarget);

figure; hold on;
title('Fitted parameter phi comparison (Camera)','fontsize',15);
plot(xticksPlot,phi_camera,'o-');
xticks(xticksPlot);
xticklabels(peaks_spd_camera_compare);
xlabel('Peak wavelength (nm)','fontsize',15);
ylabel('Fitted phi','fontsize',15);

% Add legend.
for ss = 1:length(targetCyclePerDeg)
    legendHandles{ss} = append(num2str(targetCyclePerDeg{ss}),' cpd');
end
legend(legendHandles,'fontsize',12,'location','northeastoutside');

%% Transverse Chromatic Aberration (TCA) - (SACCSFA).
%
% 1) Plot raw intensity profiles.
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Raw intensity profile over the channels (SACCSFA)');

% Make a loop to plot.
for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannelsTarget
        plot(IP_SACCSFA{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_SACCSFA(cc)),' nm');
        
        % Extract the fitted parameter, phi, for all channels and spatial
        % frequencies.
        idxParamPhi = 3;
        phi_SACCSFA(cc,ss) = params_SACCSFA{cc,ss}(idxParamPhi);
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% 2) Plot the sine fitted graphs (SACCSFA).
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Fitted intensity profile over the channels (SACCSFA)');

for ss = 1:nSFs
    subplot(5,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannelsTarget
        plot(fittedSignal_SACCSFA{cc,ss});
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% 3) Plot the comparison of the parameter phi over the channels.
figure; hold on;
title('Fitted parameter phi comparison (SACCSFA)','fontsize',15);
plot(xticksPlot,phi_SACCSFA,'o-');
xticks(xticksPlot);
xticklabels(peaks_spd_SACCSFA);
xlabel('Peak wavelength (nm)','fontsize',15);
ylabel('Fitted phi','fontsize',15);

% Add legend.
for ss = 1:length(targetCyclePerDeg)
    legendHandles{ss} = append(num2str(targetCyclePerDeg{ss}),' cpd');
end
legend(legendHandles,'fontsize',12,'location','northeastoutside');

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
    f0_lb = 29.5;
    f0_ub = 31;
    f0Range = linspace(f0_lb,f0_ub,nFits);
    
    % Set the wave to fit.
    SF = 5;
    originalSignals = IP_SACCSFA(:,SF);
    
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
