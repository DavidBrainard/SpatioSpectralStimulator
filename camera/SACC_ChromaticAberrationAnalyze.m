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
%    SACC_GetcameraImageContrast, SACC_ContrastPrintedImage_combiLED

% History:
%    11/17/23   smo    - Wrote it to use the routine for all viewing
%                        media (SACCSFA, Print, RawProjector).
%    11/22/23   smo    - Cleared up a lot and now it is working by
%                        calculating the MTF of both camera and SACCSFA
%                        within this routine.
%    12/14/23   smo    - Included 1 cpd point to all MTF measurements.
%    01/17/23   smo    - Now we interpolate the camera MTF to estimate the
%                        MTF for any spatial frequency and wavelength.
%                        Also, as a final result, we interpolated the
%                        SACCSFA MTF as well.

%% Initialize.
clear; close all;

%% Set variables.
%
% Set spatial frequency levels.
targetCyclePerDeg = {1,3,6,9,12,18};
nSFs = length(targetCyclePerDeg);

% Get the contrast calculation method.
while 1
    optionContrastCalMethod = input('Which method to calculate contrasts? [1:Average, 2:Sinefit] \n');
    if ismember(optionContrastCalMethod,[1 2])
        break
    end
    disp('Choose one between 1 (Average) and 2 (Sinefit)!');
end
switch optionContrastCalMethod
    case 1
        contrastCalMethod = 'Average';
    case 2
        contrastCalMethod = 'Sinefit';
end
fprintf('\t Contrast calculation will be based on this method - (%s) \n',contrastCalMethod);

% Get the SACCSFA trombone setting.
while 1
    tromboneSetting = input('Which Trombone setting to use? [1:Emmetropic, 2:170 nm, 3:185 nm] \n');
    if ismember(tromboneSetting,[1 2 3])
        break
    end
    disp('Choose one among 1 (Emmentropic), 2 (170 nm), 3 (185 nm)!');
end
switch tromboneSetting
    case 1
        viewingMediaSACCSFA = 'SACCSFA';
    case 2
        viewingMediaSACCSFA = 'SACCSFA170';
    case 3
        viewingMediaSACCSFA = 'SACCSFA185';
end
fprintf('\t Following mode will be run - (%s) \n',viewingMediaSACCSFA);

% Set additional analysis and plotting options. Set all these off will
% speed up running this routine.
DoFourierTransform = false;
PlotIntensityProfile = false;
PlotOneIntensityProfile = false;
PlotSineFitting = false;
PlotRawImage = false;
PlotPhiParam = false;

% Figure saving option temporarily. Set it to true will save the figures
% for the report in the current directory.
SAVEFIGURES = true;
savefileDir = '~/Desktop';

%% Get the peak wavelength of the Combi-LED (camera).
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'camera','ChromaticAberration','Spectra');
testFilename = 'CombiLED_Spectra.mat';
spdData = load(fullfile(testFiledir,testFilename));

% Extract black and white measurements per each channel.
spd_camera_white = spdData.spds.white;
spd_camera_black = spdData.spds.black;

% Get peak wavelengths.
peaks_spd_camera = FindPeakSpds(spd_camera_white,'verbose',false);

% Calculate the contrasts.
%
% Load CMFs.
S = [380 2 201];
load T_xyzJuddVos
T_XYZ = T_xyzJuddVos;
T_XYZ = 683*SplineCmf(S_xyzJuddVos,T_xyzJuddVos,S);

% Get XYZ values.
XYZ_camera_white = spd_camera_white'*T_XYZ';
XYZ_camera_black = spd_camera_black'*T_XYZ';

% Calculate contrasts.
contrasts_camera_PR670 = (XYZ_camera_white(:,2) - XYZ_camera_black(:,2))./(XYZ_camera_white(:,2) + XYZ_camera_black(:,2));

% Delete the contrast value of channel 6 as the black was not properly
% measured.
contrasts_camera_PR670(6) = NaN;

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

%% 1-a) Calculate the MTF using Average method (camera).
%
% Set the viewing media for the camera MTF measurement. We used the printed
% target so the images were saved in the folder 'Print'.
viewingMedia = 'Print';

% Load all images here.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'camera','ChromaticAberration',viewingMedia);
folders = dir(testFiledir);
dates = cell(1, numel(folders));

% Regular expression pattern to match dates in the folder names
datePattern = '\d{4}-\d{2}-\d{2}';

idxFolders = [];
% Loop through each folder and extract the date
for i = 1:numel(folders)
    folderName = folders(i).name;
    
    % Use regular expression to find the date pattern in the folder name
    match = regexp(folderName, datePattern, 'match');
    
    % Check if a date pattern was found
    if ~isempty(match)
        dates{i} = match{1};
        idxFolders(end+1) = i;
    end
end

% Extract only folders with the date in the name.
folders = folders(idxFolders);

% Remove empty cells.
dates = dates(~cellfun('isempty', dates));

% Sanity check.
if ~(numel(folders) == numel(dates))
    error(fprintf('Number of the folders (%d) and date strings (%d) does not match!',...
        numel(folders),numel(dates)));
    
end

% Get the most recent date folder directory.
dateNumbers = datenum(dates, 'yyyy-mm-dd');
[recentDateNumber, idxRecentDate] = max(dateNumbers);
recentFolderName = folders(idxRecentDate).name;
recentTestFiledir = fullfile(testFiledir,recentFolderName);

% Print out which data will be loaded.
fprintf('The data of (%s) now loading was measured on (%s) \n',viewingMedia,recentFolderName);

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
nChannels_camera = length(channelOptions);
for cc = 1:nChannels_camera
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % Make a new figure if we plot the intensity profile.
    if (PlotIntensityProfile)
        figure;
        sgtitle(sprintf('%d nm (%s)',peaks_spd_camera(cc),viewingMedia),'fontsize',15);
    end
    
    % Get the images of all spatial frequency.
    for ss = 1:nSFs
        testFilenameTemp = GetMostRecentFileName(oneChannelFileDir,...
            append(num2str(targetCyclePerDeg{ss}),'cpd_crop'));
        
        % We save all images here.
        images_camera{cc,ss} = imread(testFilenameTemp);
        
        % Get intensity profile of the image.
        image_temp = images_camera{cc,ss};
        [Ypixel Xpixel] = size(image_temp);
        
        % We will use the average of the 25% / 50% / 75% positions of the cropped image.
        IP_camera_25{cc,ss} = image_temp(round(0.25*Ypixel),:);
        IP_camera_50{cc,ss} = image_temp(round(0.50*Ypixel),:);
        IP_camera_75{cc,ss} = image_temp(round(0.75*Ypixel),:);
        
        % Set min distance between adjacent peaks.
        SF = targetCyclePerDeg{ss};
        switch SF
            % 1 cpd
            case 1
                minPeakDistance = 35;
                % 3 cpd
            case 3
                minPeakDistance = 35;
                % 6 cpd
            case 6
                minPeakDistance = 20;
            otherwise
                minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        if (PlotIntensityProfile)
            subplot(nSFs,1,ss);
            title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
        end
        
        % Calculate contrasts.
        contrastsAvg_camera_25(cc,ss) = GetIPContrast(IP_camera_25{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
        contrastsAvg_camera_50(cc,ss) = GetIPContrast(IP_camera_50{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
        contrastsAvg_camera_75(cc,ss) = GetIPContrast(IP_camera_75{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
    end
end

% Calculate the mean contrasts here.
contrastsAvg_camera = (contrastsAvg_camera_25 + contrastsAvg_camera_50 + contrastsAvg_camera_75)/3;

% Show the image if you want.
if (PlotRawImage)
    figure; hold on;
    figurePosition = [0 0 1000 1000];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('Camera captured raw images (%s)',viewingMedia));
    for cc = 1:nChannels_camera
        for ss = 1:nSFs
            subplot(nChannels_camera,nSFs,ss+nSFs*(cc-1));
            imshow(images_camera{cc,ss});
            title(sprintf('%d nm / %d cpd',peaks_spd_camera(cc),targetCyclePerDeg{ss}));
        end
    end
end

%% 1-b) Calculate the MTF using Sine fitting (camera).
%
% Here we make a loop to fit all 25%, 50%, 75% vertical positions of the
% intensity profiles
IPOptions = {'25','50','75'};
nIPOptions = length(IPOptions);

for ii = 1:nIPOptions
    
    % Set which intensity profile to load.
    whichIP = IPOptions{ii};
    
    % Load intensity profile according to the vertical positions.
    switch whichIP
        case '25'
            IP_camera = IP_camera_25;
        case '50'
            IP_camera = IP_camera_50;
        case '75'
            IP_camera = IP_camera_75;
    end
    
    
    % Load the saved initial frequencies. In not, we will newly search the
    % values.
    testFilename = fullfile(recentTestFiledir,sprintf('f0Options_%s.mat',whichIP));
    
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
            
            for cc = 1:nChannels_camera
                % Search the initial frequency here.
                signalToFit = IP_camera{cc,ss};
                f0_found(cc) = FindInitialFrequencyToFitSineWave(signalToFit,'SF',cyclesPerDeg,'verbose',false);
                
                % Show the fitting progress.
                fprintf('(%s) Searching initial frequency progress (Ch: %d/%d), (SF:%d cpd) \n',viewingMedia,cc,nChannels_camera,cyclesPerDeg);
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
        for cc = 1:nChannels_camera
            
            % Update initial frequency (f0) here.
            f0 = f0OptionsTemp(cc);
            
            % Fit happens here.
            signalToFit = IP_camera{cc,ss};
            [params_camera_temp{cc,ss}, fittedSignal_camera_temp{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
            
            % Clear the initial guess of frequency for next fit.
            clear f0;
        end
        
        % Show progress.
        fprintf('(%s) Sine fitting in progress - (%d/%d) \n',viewingMedia,ss,nSFs);
    end
    
    % Plot the results if you want.
    if (PlotSineFitting)
        for ss = 1:nSFs
            % Make a new figure per each spatial frequency.
            figure;
            figurePosition = [0 0 800 800];
            set(gcf,'position',figurePosition);
            sgtitle(sprintf('%d cpd (%s)',targetCyclePerDeg{ss},viewingMedia));
            
            % Loop over the channels.
            for cc = 1:nChannels_camera
                subplot(round(nChannels_camera/2),2,cc); hold on;
                title(sprintf('%d nm',peaks_spd_camera(cc)));
                xlabel('Pixel position');
                ylabel('dRGB');
                ylim([-10 230]);
                
                % Original.
                plot(IP_camera{cc,ss},'b-');
                
                % Fitted signal.
                plot(fittedSignal_camera_temp{cc,ss},'r-');
                legend('Original','Fit');
            end
        end
    end
    
    % Calculate contrast from the sine fitted cure.
    for ss = 1:nSFs
        for cc = 1:nChannels_camera
            paramsTemp = params_camera_temp{cc,ss};
            A = paramsTemp(1);
            B = paramsTemp(4);
            contrast = A/B;
            contrastsFit_camera_temp(cc,ss) = contrast;
        end
    end
    
    % Save out the fitting results.
    switch whichIP
        case '25'
            params_camera_25       = params_camera_temp;
            fittedSignal_camera_25 = fittedSignal_camera_temp;
            contrastsFit_camera_25 = contrastsFit_camera_temp;
        case '50'
            params_camera_50       = params_camera_temp;
            fittedSignal_camera_50 = fittedSignal_camera_temp;
            contrastsFit_camera_50 = contrastsFit_camera_temp;
        case '75'
            params_camera_75       = params_camera_temp;
            fittedSignal_camera_75 = fittedSignal_camera_temp;
            contrastsFit_camera_75 = contrastsFit_camera_temp;
    end
end

% Make an average of the contrasts.
contrastsFit_camera = (contrastsFit_camera_25 + contrastsFit_camera_50 + contrastsFit_camera_75)/3;

%% 2-a) Calculate the MTF using Average method (SACCSFA).
%
% Set viewing media to load the images.
viewingMedia = viewingMediaSACCSFA;

% Load all images here.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'camera','ChromaticAberration',viewingMedia);
folders = dir(testFiledir);
dates = cell(1, numel(folders));

% Regular expression pattern to match dates in the folder names
datePattern = '\d{4}-\d{2}-\d{2}';

idxFolders = [];
% Loop through each folder and extract the date
for i = 1:numel(folders)
    folderName = folders(i).name;
    
    % Use regular expression to find the date pattern in the folder name
    match = regexp(folderName, datePattern, 'match');
    
    % Check if a date pattern was found
    if ~isempty(match)
        dates{i} = match{1};
        idxFolders(end+1) = i;
    end
end

% Extract only folders with the date in the name.
folders = folders(idxFolders);

% Remove empty cells.
dates = dates(~cellfun('isempty', dates));

% Sanity check.
if ~(numel(folders) == numel(dates))
    error(fprintf('Number of the folders (%d) and date strings (%d) does not match!',...
        numel(folders),numel(dates)));
    
end

% Get the most recent date folder directory.
dateNumbers = datenum(dates, 'yyyy-mm-dd');
[recentDateNumber, idxRecentDate] = max(dateNumbers);
recentFolderName = folders(idxRecentDate).name;
recentTestFiledir = fullfile(testFiledir,recentFolderName);

% Print out which data will be loaded.
fprintf('The data of (%s) now loading was measured on (%s) \n',viewingMedia,recentFolderName);

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
nChannels_SACCSFA = length(channelOptions);
for cc = 1:nChannels_SACCSFA
    oneChannelFileDir = fullfile(recentTestFiledir,channelOptions{cc});
    
    % We collect the channel index here.
    idxChannels_SACCSFA(cc) = str2num(cell2mat(regexp(channelOptions{cc},'\d+','match')));
    
    % Make a new figure if we plot the intensity profile.
    if (PlotIntensityProfile)
        figure;
        sgtitle(sprintf('%d nm (%s)',peaks_spd_SACCSFA(idxChannels_SACCSFA(cc)),viewingMedia),'fontsize',15);
    end
    
    % Get the images of all spatial frequency.
    for ss = 1:nSFs
        testFilenameTemp = GetMostRecentFileName(oneChannelFileDir,...
            append(num2str(targetCyclePerDeg{ss}),'cpd_crop'));
        
        % We save all images here.
        images_SACCSFA{cc,ss} = imread(testFilenameTemp);
        
        % Get intensity profile of the image.
        image_temp = images_SACCSFA{cc,ss};
        [Ypixel Xpixel] = size(image_temp);
        
        % We will use the average of the 25% / 50% / 75% positions of the cropped image.
        IP_SACCSFA_25{cc,ss} = image_temp(round(0.25*Ypixel),:);
        IP_SACCSFA_50{cc,ss} = image_temp(round(0.50*Ypixel),:);
        IP_SACCSFA_75{cc,ss} = image_temp(round(0.75*Ypixel),:);
        
        % Set min distance between adjacent peaks.
        SF = targetCyclePerDeg{ss};
        switch SF
            % 1 cpd
            case 1
                minPeakDistance = 35;
                % 3 cpd
            case 3
                minPeakDistance = 35;
                %  cpd
            case 6
                minPeakDistance = 20;
            otherwise
                minPeakDistance = 5;
        end
        
        % Make a subplot per each spatial frequency.
        if (PlotIntensityProfile)
            subplot(nSFs,1,ss);
            title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
        end
        
        % Calculate contrasts.
        contrastsAvg_SACCSFA_25(cc,ss) = GetIPContrast(IP_SACCSFA_25{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
        contrastsAvg_SACCSFA_50(cc,ss) = GetIPContrast(IP_SACCSFA_50{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
        contrastsAvg_SACCSFA_75(cc,ss) = GetIPContrast(IP_SACCSFA_75{cc,ss},'minPeakDistance',minPeakDistance,'verbose',PlotIntensityProfile);
    end
end

% Calculate the mean contrasts here.
contrastsAvg_SACCSFA = (contrastsAvg_SACCSFA_25 + contrastsAvg_SACCSFA_50 + contrastsAvg_SACCSFA_75)/3;

% Sort the contrasts in an ascending order of the channels.
peaks_spd_SACCSFA_test = peaks_spd_SACCSFA(numChannelsSorted);
[peaks_spd_SACCSFA_test I] = sort(peaks_spd_SACCSFA_test,'ascend');
contrastsAvg_SACCSFA = contrastsAvg_SACCSFA(I,:);
IP_SACCSFA_25 = IP_SACCSFA_25(I,:);
IP_SACCSFA_50 = IP_SACCSFA_50(I,:);
IP_SACCSFA_75 = IP_SACCSFA_75(I,:);
images_SACCSFA = images_SACCSFA(I,:);

% Get number of channels to compare with the camera MTF.
nChannels_test = length(peaks_spd_SACCSFA_test);

% Show the image if you want.
if (PlotRawImage)
    figure; hold on;
    figurePosition = [0 0 1000 1000];
    set(gcf,'position',figurePosition);
    sgtitle(sprintf('Camera captured raw images (%s)',viewingMedia));
    for cc = 1:nChannels_test
        for ss = 1:nSFs
            subplot(nChannels_test,nSFs,ss+nSFs*(cc-1));
            imshow(images_SACCSFA{cc,ss});
            title(sprintf('%d nm / %d cpd',peaks_spd_SACCSFA_test(cc),targetCyclePerDeg{ss}));
        end
    end
end

% Plot the intensity profile if you want.
if (PlotOneIntensityProfile)
    whichChannel = 7;
    whichSF = 2;
    
    figure;
    figureSize = [0 0 1000 300];
    set(gcf,'position',figureSize);
    
    % Raw image.
    subplot(1,2,1);
    imshow(images_SACCSFA{whichChannel,whichSF});
    title(sprintf('Raw image \n (%d cpd / %d nm) - %s',...
        targetCyclePerDeg{whichSF},peaks_spd_SACCSFA_test(whichChannel),viewingMedia),'fontsize',15);
    
    % Intensity profile.
    subplot(1,2,2); hold on;
    plot(IP_SACCSFA_25{whichChannel,whichSF}, 'r-', 'LineWidth',1);
    plot(IP_SACCSFA_50{whichChannel,whichSF}, 'g-', 'LineWidth',1);
    plot(IP_SACCSFA_75{whichChannel,whichSF}, 'b-', 'LineWidth',1);
    title('Intensity profile','fontsize',15);
    xlabel('Pixel position (horizontal)','fontsize',15);
    ylabel('dRGB','fontsize',15);
    legend('25%','50%','75%');
end

%% 2-b) Calculate the MTF using Sine fitting method (SACCSFA).
%
% Here we make a loop to fit all 25%, 50%, 75% vertical positions of the
% intensity profiles
IPOptions = {'25','50','75'};
nIPOptions = length(IPOptions);

for ii = 1:nIPOptions
    
    % Set which intensity profile to load.
    whichIP = IPOptions{ii};
    
    % Load intensity profile according to the vertical positions.
    switch whichIP
        case '25'
            IP_SACCSFA = IP_SACCSFA_25;
        case '50'
            IP_SACCSFA = IP_SACCSFA_50;
        case '75'
            IP_SACCSFA = IP_SACCSFA_75;
    end
    
    % Load the saved initial frequencies. In not, we will newly search the
    % values.
    testFilename = fullfile(recentTestFiledir,sprintf('f0Options_%s.mat',whichIP));
    
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
            
            for cc = 1:nChannels_SACCSFA
                % Search the initial frequency here.
                signalToFit = IP_SACCSFA{cc,ss};
                f0_found(cc) = FindInitialFrequencyToFitSineWave(signalToFit,'SF',cyclesPerDeg,'verbose',false);
                
                % Show the fitting progress.
                fprintf('(%s) Searching initial frequency progress (Ch: %d/%d), (SF:%d cpd) \n',viewingMedia,cc,nChannels_SACCSFA,cyclesPerDeg);
            end
            
            % Save the value in the struct.
            f0Options = setfield(f0Options,f0FieldName,f0_found);
        end
        
        % Save out the found initial frequencies. We will load it to use next
        % time.
        save(testFilename,'f0Options');
        fprintf('All initial frequency settings found successfully and saved! - (%s) \n',viewingMedia);
    end
    
    % Fit sine signal.
    for ss = 1:nSFs
        for cc = 1:nChannels_SACCSFA
            % Set initial frequency for fitting sine wave.
            cyclesPerDeg = cell2mat(targetCyclePerDeg(ss));
            
            % Update initial frequency (f0) here.
            f0FieldName = append('SF_',num2str(cyclesPerDeg),'cpd');
            f0OptionsTemp = getfield(f0Options,f0FieldName);
            f0 = f0OptionsTemp(cc);
            
            % Fit happens here.
            signalToFit = IP_SACCSFA{cc,ss};
            [params_SACCSFA_temp{cc,ss}, fittedSignal_SACCSFA_temp{cc,ss}] = FitSineWave(signalToFit,'f0',f0,'verbose',false,'FFT',DoFourierTransform);
            
            % Clear the initial guess of frequency for next fit.
            clear f0;
        end
        
        % Show progress.
        fprintf('(%s) Sine fitting in progress - (%d/%d) \n',viewingMedia,ss,nSFs);
    end
    
    % Plot the results if you want.
    if (PlotSineFitting)
        for ss = 1:nSFs
            % Make a new figure per each spatial frequency.
            figure;
            figurePosition = [0 0 800 800];
            set(gcf,'position',figurePosition);
            sgtitle(sprintf('%d cpd (%s)',targetCyclePerDeg{ss},viewingMedia));
            
            % Loop over the channels.
            for cc = 1:nChannels_SACCSFA
                subplot(round(nChannels_SACCSFA/2),2,cc); hold on;
                title(sprintf('%d nm',peaks_spd_SACCSFA_test((cc))));
                xlabel('Pixel position');
                ylabel('dRGB');
                ylim([0 220]);
                
                % Original.
                plot(IP_SACCSFA{cc,ss},'b-');
                
                % Fitted signal.
                plot(fittedSignal_SACCSFA_temp{cc,ss},'r-');
                legend('Original','Fit');
            end
        end
    end
    
    % Calculate contrast from the sine fitted curve.
    for ss = 1:nSFs
        for cc = 1:nChannels_SACCSFA
            paramsTemp = params_SACCSFA_temp{cc,ss};
            A = paramsTemp(1);
            B = paramsTemp(4);
            contrast = A/B;
            contrastsFit_SACCSFA_temp(cc,ss) = contrast;
        end
    end
    
    % Save out the fitting results.
    switch whichIP
        case '25'
            params_SACCSFA_25       = params_SACCSFA_temp;
            fittedSignal_SACCSFA_25 = fittedSignal_SACCSFA_temp;
            contrastsFit_SACCSFA_25 = contrastsFit_SACCSFA_temp;
        case '50'
            params_SACCSFA_50       = params_SACCSFA_temp;
            fittedSignal_SACCSFA_50 = fittedSignal_SACCSFA_temp;
            contrastsFit_SACCSFA_50 = contrastsFit_SACCSFA_temp;
        case '75'
            params_SACCSFA_75       = params_SACCSFA_temp;
            fittedSignal_SACCSFA_75 = fittedSignal_SACCSFA_temp;
            contrastsFit_SACCSFA_75 = contrastsFit_SACCSFA_temp;
    end
end

% Make an average of the contrasts.
contrastsFit_SACCSFA = (contrastsFit_SACCSFA_25 + contrastsFit_SACCSFA_50 + contrastsFit_SACCSFA_75)/3;

%% 3-a) Plot the raw MTF (camera).
%
% Choose which way to calculate the contrast.
switch contrastCalMethod
    case 'Average'
        contrastRaw_camera = contrastsAvg_camera;
    case 'Sinefit'
        contrastRaw_camera = contrastsFit_camera;
end

% Calculate the contrast of the square wave from the sine wave.
factorSineToSqaurewave = 1/(4/pi);
contrastRaw_camera = contrastRaw_camera.*factorSineToSqaurewave;

% Plot the raw camera MTF results.
figure; clf;
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle(sprintf('Raw camera MTF (%s)',contrastCalMethod),'fontsize', 15);

for cc = 1:nChannels_camera
    subplot(2,4,cc); hold on;
    
    % camera MTF.
    contrastRawOneChannel = contrastRaw_camera(cc,:);
    plot(cell2mat(targetCyclePerDeg),contrastRawOneChannel,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % camera MTF (average method).
    if strcmp(contrastCalMethod,'Sinefit')
        
        contrastsAvg_1cpd = contrastsAvg_camera(:,1);
        plot(cell2mat(targetCyclePerDeg(1)),contrastsAvg_1cpd(cc),...
            'ko-','markeredgecolor','k','markerfacecolor','g', 'markersize',10);
        
        % Contrasts from PR670 measurements.
        plot(cell2mat(targetCyclePerDeg(1)),contrasts_camera_PR670(cc),...
            '+','markerfacecolor','g','markeredgecolor','k','linewidth',2,'markersize',11);
    end
    
    ylim([0 1.2]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_camera(cc)), 'fontsize', 15);
    
    % Add legend.
    switch contrastCalMethod
        case 'Average'
            legend('camera (Avg)','location','northeast','fontsize',8);
        case 'Sinefit'
            if cc == 1
                legend('camera (Sine)*pi/4','camera (Avg)','camera (PR670)','location','northeast','fontsize',8);
            elseif cc == 6
                legend('camera (Sine)*pi/4','camera (Avg)','location','southeast','fontsize',8);
            else
                legend('camera (Sine)*pi/4','camera (Avg)','camera (PR670)','location','southeast','fontsize',8);
            end
    end
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Raw_camera_MTF.tiff'));
end

%% 3-b) Interpolation of the camera MTF.
%
% Here we interpolate the camera MTF to estimate the MTF for any wavelength
% and spatial frequency combinations. We want to calculate the camera MTF
% at the same wavelengths that were used for measuring the SACCSFA MTF.
% This way, we can calculate an accruate inherent SACCSFA MTF.
z = contrastRaw_camera;
[r c] = size(z);
x = repmat(peaks_spd_camera',1,c);
y = repmat(cell2mat(targetCyclePerDeg),r,1);

% Check the matrix size.
if any(size(z) ~= size(x)) || any(size(z) ~= size(y))
    error('Matrix sizes does not match!');
end

% Fitting happens here.
% Create a lowess surface fit using fit
f_cameraMTF = fit([x(:), y(:)], z(:), 'lowess');

% Create a 3D plot to compare raw data and fitted surface
figure;
figureSize = [0 0 700 700];
set(gcf,'position',figureSize);

% Set the plot type of the fitted surface.
FittedSurfacePlotType = 2;

% Raw data.
l_raw = scatter3(x(:), y(:), z(:), 'bo','sizedata',20,'markerfacecolor','b');
hold on;

% Fitted surface.
switch FittedSurfacePlotType
    case 1
        l_fit = plot(f_cameraMTF);
    case 2
        [X, Y] = meshgrid(min(x(:)):1:max(x(:)), min(y(:)):1:max(y(:)));
        Z = feval(f_cameraMTF, [X(:), Y(:)]);
        l_fit = mesh(X, Y, reshape(Z, size(X)), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'interp');
end

title('Fitted surface - Camera MTF');
zlim([0 1]);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spatial frequency (cpd)','fontsize',15);
zlabel('Contrast','fontsize',15);
legend([l_raw l_fit], 'Raw Data','Fitted Surface');

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Camera_MTF_Interpolation.tiff'));
end

% Check how well we did the interpolation. Here we plot the measured camera
% MTF and the interpolated results together. For the measured camera MTF,
% it's the compensated results, which were used to interpolate it.
figure; hold on;
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle('Interpolated camera MTF','fontsize', 15);

for cc = 1:nChannels_camera
    peakSpdTemp = peaks_spd_camera(cc);
    
    subplot(2,4,cc); hold on;
    
    % Camera MTF - measured data.
    plot(cell2mat(targetCyclePerDeg),z(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % Camera MTF - interpolated.
    nSmoothPoints = 100;
    SF_smooth = linspace(min(cell2mat(targetCyclePerDeg)),max(cell2mat(targetCyclePerDeg)),nSmoothPoints);
    peakSpd_smooth = ones(length(SF_smooth),1).*peakSpdTemp;
    contrasts_smooth = feval(f_cameraMTF,[peakSpd_smooth,SF_smooth'])';
    plot(SF_smooth,contrasts_smooth,...
        'b-','color',[0 0 1 0.3],'linewidth',6);
    
    ylim([0 1.2]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_camera(cc)), 'fontsize', 15);
    
    % Add legend.
    switch contrastCalMethod
        case 'Average'
            legend('camera (Avg)','location','northeast','fontsize',8);
        case 'Sinefit'
            if cc == 1
                legend('camera (measure)','camera (intlp)','location','northeast','fontsize',8);
            else
                legend('camera (measure)','camera (intlp)','location','southeast','fontsize',8);
            end
    end
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Camera_MTF_Interpolation_performance.tiff'));
end

%% 3-c) Plot the camera and SACCSFA MTF together.
%
% Calculate contrast of square wave from the sine fit.
contrastsFit_SACCSFA_norm = contrastsFit_SACCSFA .* factorSineToSqaurewave;

% We used two different methods to calculate contrast. Choose either one to
% plot the results. It was chosen at the very beginning of this routine.
switch contrastCalMethod
    case 'Average'
        contrasts_SACCSFA = contrastsAvg_SACCSFA;
    case 'Sinefit'
        contrasts_SACCSFA = contrastsFit_SACCSFA_norm;
end

% Plot it here.
figure; clf;
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle(sprintf('Compensated MTF: camera vs. SACCSFA (%s)',viewingMediaSACCSFA),'fontsize', 15);

for cc = 1:nChannels_test
    subplot(2,round(nChannels_test)/2,cc); hold on;
    
    % Get the current LED channel.
    peakSpdTemp = peaks_spd_SACCSFA_test(cc);
    
    % SACCSFA MTF.
    plot(cell2mat(targetCyclePerDeg),contrasts_SACCSFA(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','r','markersize',10);
    
    % Camera MTF.
    %
    % We fit camera MTF for wavelength and spatial frequency, so here we
    % estimate the interpolated camera MTF that corresponds to the SACCSFA
    % MTF.
    %
    % Load the camera MTF.
    for ss = 1:nSFs
        sfTemp = targetCyclePerDeg{ss};
        contrasts_camera_test_temp(ss) = feval(f_cameraMTF,[peakSpdTemp,sfTemp]);
    end
    
    % Compensate the camera MTF with 1 cpd contrasts so that we can unit
    % contrast at 1 cpd at all wavelengths.
    contrasts_camera_test_temp_1cpd = contrasts_camera_test_temp(1);
    contrasts_camera_test(cc,:) = contrasts_camera_test_temp./contrasts_camera_test_temp_1cpd;
    
    % Plot it.
    plot(cell2mat(targetCyclePerDeg),contrasts_camera_test(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    ylim([0 1.2]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA_test(cc)), 'fontsize', 15);
    
    % Add legend.
    if cc == 1
        legend('SACCSFA','Camera','location','northeast','fontsize',8);
    else
        legend('SACCSFA','Camera','location','southeast','fontsize',8);
    end
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Compesated_MTF.tiff'));
end

%% 3-d) Calculate the inherent compensated MTF (SACCSFA).
%
% Here we divide the SACCSFA MTF by the camera MTF.
contrasts_SACCSFA_compensated = contrasts_SACCSFA./contrasts_camera_test;

% Set the contrast within the range.
contrasts_SACCSFA_compensated_norm = contrasts_SACCSFA_compensated;
maxContrast = 1;
contrasts_SACCSFA_compensated_norm(find(contrasts_SACCSFA_compensated_norm > maxContrast)) = 1;

% Make a new figure.
figure; hold on; clf;
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle(sprintf('Inherent SACCSFA MTF (%s)',viewingMediaSACCSFA),'fontsize', 15);

% Make a loop to plot the results of each channel.
for cc = 1:nChannels_test
    subplot(2,round(nChannels_test)/2,cc); hold on;
    contrasts_SACCSFA_compensated_temp = contrasts_SACCSFA_compensated(cc,:);
    contrasts_SACCSFA_compensated_norm_temp = contrasts_SACCSFA_compensated_norm(cc,:);
    
    % Plot it.
    %
    % Inherent SACCSFA MTF.
    plot(cell2mat(targetCyclePerDeg),contrasts_SACCSFA_compensated_temp,...
        'ko','markeredgecolor','r','markersize',13);
    
    % Inherent SACCSFA MTF (cut off the contrast over 1.0).
    plot(cell2mat(targetCyclePerDeg),contrasts_SACCSFA_compensated_norm_temp,...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    
    ylim([0 1.9]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA_test(cc)), 'fontsize', 15);
    legend('SACCSFA-raw','SACCSFA-norm','location','southeast','fontsize',10);
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'SACCSFA_MTF.tiff'));
end

%% 3-e) Interpolate the SACCSFA MTF.
%
% Here, we interpolate the SACCSFA MTF in the same way we did for the
% camera MTF so that we can estimate the contrast for any spatial frequency
% and wavelength.
%
% Clear the variables if they exist.
if exist('x')
    clear x;
end
if exist('y')
    clear y;
end
if exist('z')
    clear z;
end

% Set the variables for the interpolation.
z = contrasts_SACCSFA_compensated_norm;
[r c] = size(z);
x = repmat(peaks_spd_SACCSFA_test',1,c);
y = repmat(cell2mat(targetCyclePerDeg),r,1);

% Check the matrix size.
if any(size(z) ~= size(x)) || any(size(z) ~= size(y))
    error('Matrix sizes does not match!');
end

% Fitting happens here.
% Create a lowess surface fit using fit
f_SACCSFAMTF = fit([x(:), y(:)], z(:), 'lowess');

% Create a 3D plot to compare raw data and fitted surface
figure;
figureSize = [0 0 700 700];
set(gcf,'position',figureSize);

% Set the plot type of the fitted surface.
FittedSurfacePlotType = 2;

% Raw data.
l_raw = scatter3(x(:), y(:), z(:), 'ro','sizedata',20,'markerfacecolor','r');
hold on;

% Fitted surface.
switch FittedSurfacePlotType
    case 1
        l_fit = plot(f_SACCSFAMTF);
    case 2
        [X, Y] = meshgrid(min(x(:)):1:max(x(:)), min(y(:)):1:max(y(:)));
        Z = feval(f_SACCSFAMTF, [X(:), Y(:)]);
        l_fit = mesh(X, Y, reshape(Z, size(X)), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'interp');
end

title('Fitted surface - SACCSFA MTF');
zlim([0 1]);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spatial frequency (cpd)','fontsize',15);
zlabel('Contrast','fontsize',15);
legend([l_raw l_fit], 'Raw Data','Fitted Surface');

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'SACCSFA_MTF_Interpolation.tiff'));
end

% Check how well we did the interpolation. Here we plot the measured camera
% MTF and the interpolated results together. For the measured camera MTF,
% it's the compensated results, which were used to interpolate it.
figure; hold on;
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle('Interpolated SACCSFA MTF','fontsize', 15);

for cc = 1:nChannels_test
    peakSpdTemp = peaks_spd_SACCSFA_test(cc);
    
    subplot(2,5,cc); hold on;
    
    % SACCSFA MTF - measured data.
    plot(cell2mat(targetCyclePerDeg),z(cc,:),...
        'ko-','markeredgecolor','k','markerfacecolor','r', 'markersize',10);
    
    % SACCSFA MTF - interpolated.
    nSmoothPoints = 100;
    SF_smooth = linspace(min(cell2mat(targetCyclePerDeg)),max(cell2mat(targetCyclePerDeg)),nSmoothPoints);
    peakSpd_smooth = ones(length(SF_smooth),1).*peakSpdTemp;
    contrasts_smooth = feval(f_SACCSFAMTF,[peakSpd_smooth,SF_smooth'])';
    plot(SF_smooth,contrasts_smooth,...
        'r-','color',[1 0 0 0.3],'linewidth',6);
    
    ylim([0 1.2]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA_test(cc)), 'fontsize', 15);
    
    % Add legend.
    switch contrastCalMethod
        case 'Average'
            legend('SACCSFA (Avg)','location','northeast','fontsize',8);
        case 'Sinefit'
            legend('SACCSFA (measure)','SACCSFA (intlp)','location','southeast','fontsize',8);
    end
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'SACCSFA_MTF_Interpolation_performance.tiff'));
end

%% 4-a) Transverse Chromatic Aberration (TCA) - (camera).
%
% Plot raw intensity profiles.
figure; hold on;
figurePosition = [0 0 1000 1000];
set(gcf,'position',figurePosition);
sgtitle('Raw intensity profile over the channels (camera)');
minY = -20;
maxY = 245;

% Make a loop to plot.
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannels_camera
        plot(IP_camera_50{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_camera(cc)),' nm');
        
        % Extract the fitted parameter, phi, for all channels and spatial
        % frequencies.
        idxParamPhi = 3;
        phi_camera_25(cc,ss) = params_camera_25{cc,ss}(idxParamPhi);
        phi_camera_50(cc,ss) = params_camera_50{cc,ss}(idxParamPhi);
        phi_camera_75(cc,ss) = params_camera_75{cc,ss}(idxParamPhi);
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeastoutside','fontsize',8);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% Match up the period scale.
phi_camera_25(3:5,5) = phi_camera_25(3:5,5) + 2*pi;

% Calculate the mean phi.
phi_camera = (phi_camera_25 + phi_camera_50 + phi_camera_75)/3;

% Plot the sine fitted graphs (camera).
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Fitted intensity profile over the channels (camera)');

% Loop over Spatial frequency.
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Loop over Channel.
    for cc = 1:nChannels_camera
        plot(fittedSignal_camera_50{cc,ss});
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeastoutside','fontsize',8);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% Plot the comparison of the parameter phi over the channels.
if (PlotPhiParam)
    % Define the x-ticks for the plot.
    xticksPlot = linspace(1,nChannels_camera,nChannels_camera);
    
    figure; hold on;
    title('Fitted parameter phi comparison (camera)','fontsize',15);
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
end

% Calculate the phase shift in pixel.
%
% Get the number of the pixels. All signals should have the same size of
% the frame, so we pick one from the fitted signals.
numPixels = length(fittedSignal_camera_50{1,1});

% Get the amount of phase shift in pixel domain.
for ss = 1:nSFs
    for cc = 1:nChannels_camera
        idxParamf = 2;
        f_temp_25 = params_camera_25{cc,ss}(idxParamf);
        f_temp_50 = params_camera_50{cc,ss}(idxParamf);
        f_temp_75 = params_camera_75{cc,ss}(idxParamf);
        
        % Make an average.
        f_temp = mean([f_temp_25 f_temp_50 f_temp_75]);
        
        % Get phi parameter. If it's negative, set it to positive by adding one period (2 pi).
        phi_temp = mean([phi_camera_25(cc,ss) phi_camera_50(cc,ss) phi_camera_75(cc,ss)]);
        
        % Get period and phase shift in pixel here.
        period_pixel_camera(cc,ss) = numPixels/f_temp;
        phase_pixel_camera(cc,ss) = period_pixel_camera(cc,ss) * phi_temp/(2*pi);
    end
end

% Calculate the expected period to compare.
imageSizePixel = size(images_camera{1,1});
imageSizeHorizontalPixel = imageSizePixel(2);
imageSizeHorizontalDeg = PixelToDeg(imageSizeHorizontalPixel,'dir','horizontal','verbose',false);
imageSizeHorizontalPixelOneDeg = imageSizeHorizontalPixel/imageSizeHorizontalDeg;
period_pixel_camera_expected = imageSizeHorizontalPixelOneDeg./cell2mat(targetCyclePerDeg);

% Plot the period in pixel per spatial frequency.
figure;
figureSize = [0 0 450 1000];
set(gcf,'position',figureSize);
sgtitle('Sine fitted period in pixel (camera)');
x_data = linspace(1,nChannels_camera,nChannels_camera);
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Measured period.
    plot(x_data, period_pixel_camera(:,ss),'b-o','markerfacecolor','b','markeredgecolor','k');
    % Expected period.
    plot([min(x_data) max(x_data)],ones(1,2).*period_pixel_camera_expected(ss),'b-','color',[0 0 1 0.3],'linewidth',3);
    
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_camera);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Period (pixel)','fontsize',15);
    ylim([0 period_pixel_camera_expected(ss)*2+1]);
    yticks(round([0 period_pixel_camera_expected(ss) period_pixel_camera_expected(ss)*2]));
    legend('Measure','Expected');
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Camera_period_pixel.tiff'));
end

% Plot the phase shift in pixel per spatial frequency.
%
% We will compare based on the channel that we focused with the camera.
channelFocus = 598;
idxChannelFocus = find(peaks_spd_camera == channelFocus);
phase_pixel_camera_ref = phase_pixel_camera(idxChannelFocus,:);
phaseShift_pixel_camera = phase_pixel_camera_ref - phase_pixel_camera;

% Plot the phase shift in pixel.
figure;
figureSize = [0 0 450 1000];
set(gcf,'position',figureSize);

sgtitle('Phase shift in pixel (camera)');
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Measured phase shift.
    plot(x_data, phaseShift_pixel_camera(:,ss),'b-o','markerfacecolor','b','markeredgecolor','k');
    % No difference line.
    plot([min(x_data) max(x_data)], zeros(1,2),'b-','color',[0 0 1 0.3],'linewidth',3);
    
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticklabels(peaks_spd_camera);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Shift (pixel)','fontsize',15);
    ylim([-5 5]);
    legend('Measure','No difference');
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'Camera_PhaseShift_pixel.tiff'));
end

% Calculate phase shift in degrees and the phase shift in pixel on the DMD
% of the SACCSFA system. We use the function PixelToDeg which will
% calculate both.
[phaseShift_deg_camera phaseShift_pixel_camera_DMD] = PixelToDeg(phaseShift_pixel_camera,'verbose',false);

% Get the fitted spatial frequency.
fittedSF_camera = 1./PixelToDeg(period_pixel_camera,'verbose',false);

%% 4-b) Transverse Chromatic Aberration (TCA) - (SACCSFA).
%
% Plot raw intensity profiles.
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Raw intensity profile over the channels (SACCSFA)');

% Make a loop to plot.
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannels_test
        plot(IP_SACCSFA_50{cc,ss});
        
        % Generate texts for the legend for each graph.
        legendHandles{cc} = append(num2str(peaks_spd_SACCSFA_test(cc)),' nm');
        
        % Extract the fitted parameter, phi, for all channels and spatial
        % frequencies.
        idxParamPhi = 3;
        phi_SACCSFA_25(cc,ss) = params_SACCSFA_25{cc,ss}(idxParamPhi);
        phi_SACCSFA_50(cc,ss) = params_SACCSFA_50{cc,ss}(idxParamPhi);
        phi_SACCSFA_75(cc,ss) = params_SACCSFA_75{cc,ss}(idxParamPhi);
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeastoutside','fontsize',8);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% Correct phi to calculate the phase shift correct. For now, we
% manually correct it, but maybe we want to do this part more
% elaborately later on.
switch viewingMediaSACCSFA
    case 'SACCSFA'
        phi_SACCSFA_25(10,3) = phi_SACCSFA_25(10,3) - 2*pi;
    case 'SACCSFA170'
        phi_SACCSFA_25(1,4) = phi_SACCSFA_25(1,4) - 2*pi;
        phi_SACCSFA_50(1,4) = phi_SACCSFA_50(1,4) - 2*pi;
        phi_SACCSFA_75(1,4) = phi_SACCSFA_75(1,4) - 2*pi;
        
        phi_SACCSFA_25(2,4) = phi_SACCSFA_25(2,4) - 2*pi;
        phi_SACCSFA_50(2,4) = phi_SACCSFA_50(2,4) - 2*pi;
        phi_SACCSFA_75(2,4) = phi_SACCSFA_75(2,4) - 2*pi;
        
        phi_SACCSFA_75(9,4) = phi_SACCSFA_75(9,4) - 2*pi;
    case 'SACCSFA185'
        phi_SACCSFA_50(1,5) = phi_SACCSFA_50(1,5) - 2*pi;
        phi_SACCSFA_75(1,5) = phi_SACCSFA_75(1,5) - 2*pi;
        
        phi_SACCSFA_50(2,5) = phi_SACCSFA_50(2,5) - 2*pi;
end

% Calculate the mean phi here.
phi_SACCSFA = (phi_SACCSFA_25 + phi_SACCSFA_50 + phi_SACCSFA_75)/3;

% Plot the sine fitted graphs (SACCSFA).
figure; hold on;
set(gcf,'position',figurePosition);
sgtitle('Fitted intensity profile over the channels (SACCSFA)');

for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Channel.
    for cc = 1:nChannels_test
        plot(fittedSignal_SACCSFA_50{cc,ss});
    end
    
    % Set each graph in the same format.
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    legend(legendHandles,'fontsize',11,'location','southeastoutside','fontsize',8);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('dRGB','fontsize',12);
    ylim([minY maxY]);
end

% Plot the comparison of the parameter phi over the channels.
if (PlotPhiParam)
    % Define the x-ticks for the plot.
    xticksPlot = linspace(1,nChannels_test,nChannels_test);
    
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
end

% Calculate the phase shift in pixel.
%
% Get the number of the pixels. All signals should have the same size
% of the frame, so we pick one from the fitted signals.
numPixels = length(fittedSignal_SACCSFA_50{1,1});

% Get the amount of phase shift in pixel domain.
for ss = 1:nSFs
    % Get spatial frequency.
    SF = targetCyclePerDeg{ss};
    
    for cc = 1:nChannels_test
        
        % Get frequency.
        idxParamf = 2;
        f_temp_25 = params_SACCSFA_25{cc,ss}(idxParamf);
        f_temp_50 = params_SACCSFA_50{cc,ss}(idxParamf);
        f_temp_75 = params_SACCSFA_75{cc,ss}(idxParamf);
        
        % Make an average.
        f_temp = mean([f_temp_25 f_temp_50 f_temp_75]);
        
        % Get phi parameter. If it's negative, set it to positive by adding one period (2 pi).
        phi_temp = mean([phi_SACCSFA_25(cc,ss) phi_SACCSFA_50(cc,ss) phi_SACCSFA_75(cc,ss)]);
        
        % Get period and phase shift in pixel here.
        period_pixel_SACCSFA(cc,ss) = numPixels/f_temp;
        
        % Calculate the phase shift in pixel here.
        phase_pixel_SACCSFA(cc,ss) = period_pixel_SACCSFA(cc,ss) * phi_temp/(2*pi);
    end
end

% Calculate the expected period to compare.
imageSizePixel = size(images_SACCSFA{1,1});
imageSizeHorizontalPixel = imageSizePixel(2);
imageSizeHorizontalDeg = PixelToDeg(imageSizeHorizontalPixel,'dir','horizontal','verbose',false);
imageSizeHorizontalPixelOneDeg = imageSizeHorizontalPixel/imageSizeHorizontalDeg;
period_pixel_SACCSFA_expected = imageSizeHorizontalPixelOneDeg./cell2mat(targetCyclePerDeg);

% Plot the period in pixel per channel.
figure;
figureSize = [0 0 450 1000];
set(gcf,'position',figureSize);
sgtitle(sprintf('Sine fitted period in pixel (%s)',viewingMediaSACCSFA));

x_data = linspace(1,nChannels_test,nChannels_test);
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Measured period.
    plot(x_data, period_pixel_SACCSFA(:,ss),'r-o','markerfacecolor','r','markeredgecolor','k');
    % Expected period.
    plot([min(x_data) max(x_data)], ones(1,2).*period_pixel_SACCSFA_expected(ss),'r-','color',[1 0 0 0.3],'linewidth',3);
    
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticks(x_data);
    xticklabels(peaks_spd_SACCSFA_test);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Period (pixel)','fontsize',15);
    ylim([0 period_pixel_SACCSFA_expected(ss)*2+1]);
    yticks(round([0 period_pixel_SACCSFA_expected(ss) period_pixel_SACCSFA_expected(ss)*2]));
    legend('Measure','Expected');
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'SACCSFA_period_pixel.tiff'));
end

% Plot the phase shift in pixel per spatial frequency.
%
% We will compare based on the channel that we focused with the camera.
channelFocus = 592;
idxChannelFocus = find(peaks_spd_SACCSFA_test == channelFocus);
phase_pixel_SACCSFA_ref = phase_pixel_SACCSFA(idxChannelFocus,:);
phaseShift_pixel_SACCSFA = phase_pixel_SACCSFA_ref - phase_pixel_SACCSFA;

% Plot happens here.
figure;
figureSize = [0 0 450 1000];
set(gcf,'position',figureSize);

sgtitle(sprintf('Phase shift in pixel (%s)',viewingMediaSACCSFA));
for ss = 1:nSFs
    subplot(nSFs,1,ss); hold on;
    
    % Measured phase shift.
    plot(x_data, phaseShift_pixel_SACCSFA(:,ss),'r-o','markerfacecolor','r','markeredgecolor','k');
    % No difference line.
    plot([min(x_data) max(x_data)], zeros(1,2),'r-','color',[1 0 0 0.3],'linewidth',3);
    
    title(sprintf('%d cpd',targetCyclePerDeg{ss}),'fontsize',15);
    xticks(x_data);
    xticklabels(peaks_spd_SACCSFA_test);
    xlabel('Peak wavelength (nm)','fontsize',15);
    ylabel('Shift (pixel)','fontsize',15);
    ylim([-5 5]);
    legend('Measure','No difference');
end

% Save the image on the Desktop if you want.
if (SAVEFIGURES)
    saveas(gcf,fullfile(savefileDir,'SACCSFA_PhaseShift_pixel.tiff'));
end

% Calculate phase shift in degrees and the phase shift in pixel on the DMD
% of the SACCSFA system. We use the function PixelToDeg which will
% calculate both.
[phaseShift_deg_SACCSFA phaseShift_pixel_SACCSFA_DMD] = PixelToDeg(phaseShift_pixel_SACCSFA,'verbose',false);

% Get the fitted spatial frequency.
fittedSF_SACCSFA = 1./PixelToDeg(period_pixel_SACCSFA,'verbose',false);

%% Plot the channels that we used in this study.
PLOTSPECTRUM = false;

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
    
    % camera.
    figure; hold on;
    p3 = plot(wls,spd_camera_white,'k-');
    p4 = plot(wls,spd_camera_white,'-','linewidth',5,'color',[0 1 0 0.3]);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    ylim([0 max(spd_camera_white,[],'all')*1.01]);
    legend([p3(1) p4(1)],'All channels (N=8)','Tested channel (N=8)','fontsize',12,'location','northeast');
    title('camera (Combi-LED)','fontsize',15);
end
