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

%% Get the peak wavelength of the Combi-LED (Camera).
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration','Spectra');
testFilename = 'CombiLED_Spectra.mat';
spdData = load(fullfile(testFiledir,testFilename));

% Extract black and white measurements per each channel.
spd_camera = spdData.spds.white;
spd_camera_black = spdData.spds.black;

% Get peak wavelengths.
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

% Sort the contrasts in an ascending order of the channels.
peaks_spd_SACCSFA_test = peaks_spd_SACCSFA(numChannelsSorted);
[peaks_spd_SACCSFA_test I] = sort(peaks_spd_SACCSFA_test,'ascend');
meanContrasts_SACCSFA = meanContrasts_SACCSFA(I,:);
IP_SACCSFA = IP_SACCSFA(I,:);

% Get number of channels to compare with the camera MTF.
nChannelsTest = length(peaks_spd_SACCSFA_test);

%% Fit sine function to the signal (SACCSFA).
%
% Load the saved initial frequencies. In not, we will newly search the
% values.
testFilename = fullfile(recentTestFiledir,'f0Options.mat');

% Load the file here.
if isfile(testFilename)
    fprintf('Pre-saved f0 options file found! We will load it for sine fitting - (%s) \n',viewingMedia);
    f0Options = load(testFilename);
    f0Options = f0Options.f0Options;
else
    % If not, we will search the initial frequencies to fit sine to the
    % intensity profiles.
    fprintf('We will start searching for initial frequency (f0) for sine fitting - (%s) \n',viewingMedia);
    f0Options = struct();
    for ss = 1:nSFs
        % Set spatial frequency.
        cyclesPerDeg = cell2mat(targetCyclePerDeg(ss));
        
        % Set field name in the struct.
        f0FieldName = append('SF_',num2str(cyclesPerDeg),'cpd');
        
        for cc = 1:nChannels
            % Search the initial frequency here.
            signalToFit = IP_SACCSFA{cc,ss};
            f0_found(cc) = FindInitialFrequencyToFitSineWave(signalToFit,'SF',cyclesPerDeg,'verbose',false);
            
            % Show the fitting progress.
            fprintf('(%s) Searching initial frequency progress (Ch: %d/%d), (SF:%d cpd) \n',viewingMedia,cc,nChannels,cyclesPerDeg);
        end
        
        % Save the value in the struct.
        f0Options = setfield(f0Options,f0FieldName,f0_found);
    end
    
    % Save out the found initial frequencies. We will load it to use next
    % time.
    save(testFilename,'f0Options');
end

% Fit sine signal.
for ss = 1:nSFs
    for cc = 1:nChannels
        % Set initial frequency for fitting sine wave.
        cyclesPerDeg = cell2mat(targetCyclePerDeg(ss));
        
        % Update initial frequency (f0) here.
        f0FieldName = append('SF_',num2str(cyclesPerDeg),'cpd');
        f0OptionsTemp = getfield(f0Options,f0FieldName);
        f0 = f0OptionsTemp(cc);
        
        % Fit happens here.
        signalToFit = IP_SACCSFA{cc,ss};
        [params_SACCSFA{cc,ss}, fittedSignal_SACCSFA{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
    
    % Show progress.
    fprintf('(%s) Sine fitting in progress - (%d/%d) \n',viewingMedia,ss,nSFs);
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
        title(sprintf('%d nm',peaks_spd_SACCSFA_test((cc))));
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
        contrastsFit_SACCSFA(cc,ss) = contrast;
    end
end

%% Get contrasts of 1 cpd image for compensation (SACCSFA).
figure; hold on;
sgtitle(sprintf('1 cpd (%s)',viewingMedia));
for cc = 1:nChannelsTest
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    testFilename = '1cpd_crop';
    testFilename = GetMostRecentFileName(oneChannelFileDir,testFilename);
    
    % We save an image here.
    image = imread(testFilename);
    
    % Set min distance between adjacent peaks.
    minPeakDistance = 30;
    
    % Calculate contrasts.
    if (plotIntensityProfile)
        figure;
        title(sprintf('%d nm (%s)',peaks_spd_SACCSFA_test(cc),viewingMedia),'fontsize',15);
        subtitle('1 cpd','fontsize',13);
    end
    [contrastsTemp, IP_SACCSFA_1cpd{cc}] = GetImgContrast(image,'minPeakDistance',minPeakDistance,'verbose',plotIntensityProfile);
    meanContrasts_SACCSFA_1cpd(cc) = mean(contrastsTemp);
    
    % Now fit sine curve to the 1 cpd to calculate contrast.
    %
    % Set initial f0 options for sine fitting.
    f0Options = [3.275862 3.275862 3.275862 3.275862 3.275862 3.379310 3.275862 3.379310 3.275862 3.379310];
    f0 = f0Options(cc);
    
    % Fit happens here.
    signalToFit = IP_SACCSFA_1cpd{cc};
    [params_SACCSFA_1cpd{cc}, fittedSignal_SACCSFA_1cpd{cc}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
    
    % Clear the initial guess of frequency for next fit.
    clear f0;
    
    % Plot the results.
    subplot(round(nChannels/2),2,cc); hold on;
    title(sprintf('%d nm',peaks_spd_SACCSFA_test(cc)));
    xlabel('Pixel position');
    ylabel('dRGB');
    ylim([-20 240]);
    
    % Original.
    plot(IP_SACCSFA_1cpd{cc},'b-');
    
    % Fitted signal.
    plot(fittedSignal_SACCSFA_1cpd{cc},'r-');
    legend('Origianl','Fit');
    
    % Calculate contrast.
    paramsTemp = params_SACCSFA_1cpd{cc};
    A = paramsTemp(1);
    B = paramsTemp(4);
    contrast = A/B;
    contrastsFit_SACCSFA_1cpd(cc) = contrast;
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
end

%% Fit sine function to the signal (Camera).
%
% Load the saved initial frequencies. In not, we will newly search the
% values.
testFilename = fullfile(recentTestFiledir,'f0Options.mat');

% Load the file here.
if isfile(testFilename)
    fprintf('Pre-saved f0 options file found! We will load it for sine fitting - (%s) \n',viewingMedia);
    f0Options = load(testFilename);
    f0Options = f0Options.f0Options;
else
    % If not, we will search the initial frequencies to fit sine to the
    % intensity profiles.
    fprintf('We will start searching for initial frequency (f0) for sine fitting - (%s) \n',viewingMedia);
    f0Options = struct();
    for ss = 1:nSFs
        % Set spatial frequency.
        cyclesPerDeg = cell2mat(targetCyclePerDeg(ss));
        
        % Set field name in the struct.
        f0FieldName = append('SF_',num2str(cyclesPerDeg),'cpd');
        
        for cc = 1:nChannels
            % Search the initial frequency here.
            signalToFit = IP_SACCSFA{cc,ss};
            f0_found(cc) = FindInitialFrequencyToFitSineWave(signalToFit,'SF',cyclesPerDeg,'verbose',false);
            
            % Show the fitting progress.
            fprintf('(%s) Searching initial frequency progress (Ch: %d/%d), (SF:%d cpd) \n',viewingMedia,cc,nChannels,cyclesPerDeg);
        end
        
        % Save the value in the struct.
        f0Options = setfield(f0Options,f0FieldName,f0_found);
    end
    
    % Save out the found initial frequencies. We will load it to use next
    % time.
    save(testFilename,'f0Options');
end

% Fit sine signal here.
%
% Loop over the spatial frequency.
for ss = 1:nSFs
    % Set initial frequency for fitting sine wave.
    cyclesPerDeg = cell2mat(targetCyclePerDeg(ss));
    f0FieldName = append('SF_',num2str(cyclesPerDeg),'cpd');
    f0OptionsTemp = getfield(f0Options,f0FieldName);
    
    % Loop over the channels.
    for cc = 1:nChannels
        
        % Update initial frequency (f0) here.
        f0 = f0OptionsTemp(cc);
        
        % Fit happens here.
        signalToFit = IP_camera{cc,ss};
        [params_camera{cc,ss}, fittedSignal_camera{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
        
        % Clear the initial guess of frequency for next fit.
        clear f0;
    end
    
    % Show progress.
    fprintf('(%s) Sine fitting in progress - (%d/%d) \n',viewingMedia,ss,nSFs);
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
        contrastsFit_camera(cc,ss) = contrast;
    end
end

%% Get contrasts of 1 cpd image for compensation (Camera).
figure; hold on;
for cc = 1:nChannels
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    testFilename = '1cpd_crop';
    testFilename = GetMostRecentFileName(oneChannelFileDir,testFilename);
    
    % We save an image here.
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
    %
    % Set initial f0 options for sine fitting.
    f0Options = [2.551724 2.551724 2.551724 2.551724 2.551724 2.551724 2.551724 1.931034];
    f0 = f0Options(cc);
    
    % Fit happens here.
    signalToFit = IP_camera_1cpd{cc};
    [params_camera_1cpd{cc}, fittedSignal_camera_1cpd{cc}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
    
    % Clear the initial guess of frequency for next fit.
    clear f0;
    
    % Plot the results.
    subplot(round(nChannels/2),2,cc); hold on;
    title(sprintf('%d nm',peaks_spd_camera(cc)));
    xlabel('Pixel position');
    ylabel('dRGB');
    ylim([-20 240]);
    
    % Original.
    plot(IP_camera_1cpd{cc},'b-');
    
    % Fitted signal.
    plot(fittedSignal_camera_1cpd{cc},'r-');
    legend('Origianl','Fit');
    
    % Calculate contrast.
    paramsTemp = params_camera_1cpd{cc};
    A = paramsTemp(1);
    B = paramsTemp(4);
    contrast = A/B;
    contrastsFit_camera_1cpd(cc) = contrast;
end

%% 2) Plot the raw MTF and compensate it (Camera).
%
% Choose which way to calculate the contrast.
switch contrastCalMethod
    case 'MeanIntensityProfile'
        contrastRaw_camera = meanContrasts_camera;
    case 'Sinefit'
        contrastRaw_camera = contrastsFit_camera;
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
    
    legend('Camera (Raw)','location','southeast','fontsize',11);
    ylim([0 1.15]);
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
meanContrasts_cameraNorm = meanContrasts_camera./meanContrasts_camera_1cpd';
contrastsFit_cameraNorm = contrastsFit_camera./contrastsFit_camera_1cpd';

meanContrasts_SACCSFANorm = meanContrasts_SACCSFA./meanContrasts_SACCSFA_1cpd';
contrastsFit_SACCSFANorm = contrastsFit_SACCSFA./contrastsFit_SACCSFA_1cpd';

%% 3) Calculate the compensated MTF (SACCSFA).
%
% We used two different methods to calculate contrast. Choose either one to
% plot the results. It was chosen at the very beginning of this routine.
switch contrastCalMethod
    case 'MeanIntensityProfile'
        contrast_camera = meanContrasts_cameraNorm;
        contrast_SACCSFA = meanContrasts_SACCSFANorm;
    case 'Sinefit'
        contrast_camera = contrastsFit_cameraNorm;
        contrast_SACCSFA = contrastsFit_SACCSFANorm;
end

% Here we choose which channel of the combi-LED to compare to each channel
% of SACCSFA. We will choose the one with the closest peak wavelength
% between combi-LED and SACCSFA to compare.
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
sgtitle('Compensated MTF: Camera vs. SACCSFA', 'fontsize', 15);

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
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % Save out camera MTF to test.
    %
    % We will use this to compensate the SACCSFA MTF later on.
    contrast_camera_test(cc,:) = meanContrastsOneChannel;
    
    % Plot stuffs.
    legend(sprintf('SACCSFA (%d nm)',peaks_spd_SACCSFA_test(cc)),...
        sprintf('Camera (%d nm)',peaks_spd_camera_test(cc)),'location','southeast','fontsize',8);
    ylim([0 1.15]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA_test(cc)), 'fontsize', 15);
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
sgtitle('Compensated SACCSFA MTF (SACCSFA MTF/Camera MTF)', 'fontsize', 15);

% Make a loop to plot the results of each channel.
for cc = 1:nChannelsTest
    subplot(2,round(nChannelsTest)/2,cc); hold on;
    meanContrastsOneChannel = contrast_SACCSFA_compensated(cc,:);
    plot(cell2mat(targetCyclePerDeg),meanContrastsOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    ylim([0 1.15]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA_test(cc)), 'fontsize', 15);
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
maxY = 230;

% Make a loop to plot.
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannels
        plot(IP_camera{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_camera(cc)),' nm');
        
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
    subplot(nSFs,1,ss); hold on;
    
    % Loop over Channel.
    for cc = 1:nChannels
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
xticksPlot = linspace(1,nChannels,nChannels);

figure; hold on;
title('Fitted parameter phi comparison (Camera)','fontsize',15);
plot(xticksPlot,phi_camera,'o-');
xticks(xticksPlot);
xticklabels(peaks_spd_camera);
xlabel('Peak wavelength (nm)','fontsize',15);
ylabel('Fitted phi','fontsize',15);

% Add legend.
clear legendHandles;
for ss = 1:length(targetCyclePerDeg)
    legendHandles{ss} = append(num2str(targetCyclePerDeg{ss}),' cpd');
end
legend(legendHandles,'fontsize',12,'location','northeastoutside');

% Calculate the phase shift in pixel.
%
% Get the number of the pixels. All signals should have the same size of
% the frame, so we pick one from the fitted signals.
numPixels = length(fittedSignal_camera{1,1});

% Get the amount of phase shift in pixel domain.
for ss = 1:nSFs
    for cc = 1:nChannels
        params_temp = params_camera{cc,ss};
        fitted_f_temp = params_temp(2);
        
        % Get phi parameter. If it's negative, set it to positive by adding one period (2 pi).
        fitted_phi_temp = params_temp(3);
        phi_camera(cc,ss) = fitted_phi_temp;
        
        % Get period and phase shift in pixel here.
        onePeriod_pixel_camera(cc,ss) = numPixels/fitted_f_temp;
        phaseShift_pixel_camera(cc,ss) = onePeriod_pixel_camera(cc,ss) * fitted_phi_temp/(2*pi);
    end
end

% Plot the period in pixel per spatial frequency..
figure;
figureSize = [0 0 600 800];
set(gcf,'position',figureSize);
sgtitle('Sine fitted period in pixel (Camera)');
x_data = linspace(1,nChannels,nChannels);
for ss = 1:nSFs
    subplot(nSFs,1,ss);
    plot(x_data, onePeriod_pixel_camera(:,ss),'b-o','markerfacecolor','b','markeredgecolor','k');
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_camera);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Period (pixel)','fontsize',15);
    ylim([0 1.3*max(onePeriod_pixel_camera,[],'all')]);
end

% Plot the phase shift in pixel per spatial frequency.
%
% We will compare based on the channel that we focused with the camera.
channelFocus = 598;
idxChannelFocus = find(peaks_spd_camera == channelFocus);
phaseShift_pixel_camera_ref = phaseShift_pixel_camera(idxChannelFocus,:);
phaseShift_pixel_camera_diff = abs(round(phaseShift_pixel_camera - phaseShift_pixel_camera_ref,1));

% Plot the phase shift in pixel.
figure;
figureSize = [0 0 600 800];
set(gcf,'position',figureSize);

sgtitle('Phase shift in pixel (Camera)');
for ss = 1:nSFs
    subplot(nSFs,1,ss);
    plot(x_data, phaseShift_pixel_camera_diff(:,ss),'b-o','markerfacecolor','b','markeredgecolor','k');
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_camera);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Shift (pixel)','fontsize',15);
    ylim([0 4]);
end

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
    for cc = 1:nChannelsTest
        plot(IP_SACCSFA{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_SACCSFA_test(cc)),' nm');
        
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
    for cc = 1:nChannelsTest
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
%
% Define the x-ticks for the plot.
xticksPlot = linspace(1,nChannelsTest,nChannelsTest);

figure; hold on;
title('Fitted parameter phi comparison (SACCSFA)','fontsize',15);
plot(xticksPlot,phi_SACCSFA,'o-');
xticks(xticksPlot);
xticklabels(peaks_spd_SACCSFA_test);
xlabel('Peak wavelength (nm)','fontsize',15);
ylabel('Fitted phi','fontsize',15);

% Add legend.
clear legendHandles;
for ss = 1:length(targetCyclePerDeg)
    legendHandles{ss} = append(num2str(targetCyclePerDeg{ss}),' cpd');
end
legend(legendHandles,'fontsize',12,'location','northeastoutside');

% Calculate the phase shift in pixel.
%
% Get the number of the pixels. All signals should have the same size of
% the frame, so we pick one from the fitted signals.
numPixels = length(fittedSignal_SACCSFA{1,1});

% Get the amount of phase shift in pixel domain.
for ss = 1:nSFs
    for cc = 1:nChannelsTest
        params_temp = params_SACCSFA{cc,ss};
        fitted_f_temp = params_temp(2);
        
        % Get phi parameter. If it's negative, set it to positive by adding one period (2 pi).
        fitted_phi_temp = params_temp(3);
        phi_SACCSFA(cc,ss) = fitted_phi_temp;
        
        
        
        
        % THIS IS TEMP SOLUTION WHICH WILL BE FIXED (TEMPSOLUTION).
        if and(cc==1,ss==2)
            phi_SACCSFA(cc,ss)=  phi_SACCSFA(cc,ss) - 2*pi;
        end
        
        
        
        
        % Get period and phase shift in pixel here.
        onePeriod_pixel_SACCSFA(cc,ss) = numPixels/fitted_f_temp;
        
        % Calculate the phase shift in pixel here.
        phaseShift_pixel_SACCSFA(cc,ss) = onePeriod_pixel_SACCSFA(cc,ss) * phi_SACCSFA(cc,ss)/(2*pi);
    end
end

% Plot the period in pixel per channel.
figure;
figureSize = [0 0 600 800];
set(gcf,'position',figureSize);

sgtitle('Sine fitted period in pixel (SACCSFA)');
x_data = linspace(1,nChannelsTest,nChannelsTest);
for ss = 1:nSFs
    subplot(nSFs,1,ss);
    plot(x_data, onePeriod_pixel_SACCSFA(:,ss),'r-o','markerfacecolor','r','markeredgecolor','k');
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_SACCSFA_test);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Period (pixel)','fontsize',15);
    ylim([0 1.3*max(onePeriod_pixel_SACCSFA,[],'all')]);
end

% Plot the phase shift in pixel per spatial frequency.
%
% We will compare based on the channel that we focused with the camera.
channelFocus = 592;
idxChannelFocus = find(peaks_spd_SACCSFA_test == channelFocus);
phaseShift_pixel_SACCSFA_ref = phaseShift_pixel_SACCSFA(idxChannelFocus,:);
phaseShift_pixel_SACCSFA_diff = abs(round(phaseShift_pixel_SACCSFA - phaseShift_pixel_SACCSFA_ref,1));

% Plot the phase shift in pixel.
figure;
figureSize = [0 0 600 800];
set(gcf,'position',figureSize);

sgtitle('Phase shift in pixel (SACCSFA)');
for ss = 1:nSFs
    subplot(nSFs,1,ss);
    plot(x_data, phaseShift_pixel_SACCSFA_diff(:,ss),'r-o','markerfacecolor','r','markeredgecolor','k');
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_SACCSFA_test);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Shift (pixel)','fontsize',15);
    ylim([0 4]);
end

%% Plot the channels that we used in this study.
PLOTSPECTRUM = true;

if (PLOTSPECTRUM)
    % SACCSFA.
    figure; hold on;
    S = recentCalData.rawData.S;
    wls = SToWls(S);
    p1 = plot(wls,spd_SACCSFA,'k-');
    p2 = plot(wls,spd_SACCSFA(:,idxChannels_SACCSFA),'-','linewidth',5,'color',[1 0 0 0.3]);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    ylim([0 max(spd_SACCSFA,[],'all')*1.01]);
    legend([p1(1) p2(1)],'All channels (N=16)','Tested channel (N=10)','fontsize',12,'location','northwest');
    title('SACCSFA','fontsize',15);
    
    % Camera.
    figure; hold on;
    p3 = plot(wls,spd_camera,'k-');
    p4 = plot(wls,spd_camera,'-','linewidth',5,'color',[0 1 0 0.3]);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    ylim([0 max(spd_camera,[],'all')*1.01]);
    legend([p3(1) p4(1)],'All channels (N=8)','Tested channel (N=8)','fontsize',12,'location','northeast');
    title('Camera (Combi-LED)','fontsize',15);
    
    % Spectrum comparison: SACCSFA vs. Camera
    figure; hold on;
    for tt = 1:nChannelsTest
        subplot(2,round(nChannelsTest/2),tt); hold on;
        
        xlabel('Wavelength (nm)','fontsize',15);
        ylabel('Spectral power','fontsize',15);
        ylim([0 max(spd_camera,[],'all')*1.01]);
        
        title('Camera (Combi-LED)','fontsize',15);
    end
end
