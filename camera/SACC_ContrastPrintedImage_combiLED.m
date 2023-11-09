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
measureDate = '0905';
focusedImage = false;

%% Plot spectra.
PLOTSPECTRA = false;
measureDateSPD = '0829';

testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDateSPD);
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

% Plot it - Raw.
figure; clf;
plot(wls,spd,'linewidth',1);
xlabel('Wavelength (nm)','fontsize',15);
ylabel('Spectral power','fontsize',15);
xticks([380:80:780]);
for ll = 1:size(spd,2)
    legendHandles{ll} = append(num2str(peaks_spd(ll)),' nm');
end
legend(legendHandles,'fontsize',15);

% Plot it - Normalized.
if (PLOTSPECTRA)
    figure; clf;
    plot(wls,spd./max(spd),'linewidth',1);
    xlabel('Wavelength (nm)','fontsize',15);
    ylabel('Spectral power','fontsize',15);
    xticks([380:80:780]);
    legend(legendHandles,'fontsize',15);
end

%% Load SACCSFA MTF results to compare.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDate);
testFilename = 'MTF_SACCSFA.mat';
data_SACCSFA = load(fullfile(testFiledir,testFilename));
meanContrasts_SACCSFA = data_SACCSFA.meanContrasts_all;

%% Load all images here.
for cc = 1:8
    targetChOptions{cc} = append('Ch',num2str(cc));
end
nTargetChs = length(targetChOptions);

% Make a loop here.
for cc = 1:nTargetChs
    targetChTemp = targetChOptions{cc};
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Camera','PrintedImageContrast',measureDate,targetChTemp);
    else
        error('Cannot find data file list!');
    end
    
    nSFs = length(targetCyclePerDeg);
    for dd = 1:nSFs
        % Get the file name of the images.
        if (focusedImage)
            % Load focused image if there is any available. If not, this
            % line will be skipped.
            try
                testFilenameTemp = GetMostRecentFileName(testFiledir,...
                    append(num2str(targetCyclePerDeg{dd}),'cpd_focused_crop'));
            catch
                % Alarm message.
                disp('There is no such file, so regular image will be loaded');
                
                % Load the regular image measured with fixed focus (infinity).
                testFilenameTemp = GetMostRecentFileName(testFiledir,...
                    append(num2str(targetCyclePerDeg{dd}),'cpd_crop'));
            end
        else
            % Load the regular image measured with fixed focus (infinity).
            testFilenameTemp = GetMostRecentFileName(testFiledir,...
                append(num2str(targetCyclePerDeg{dd}),'cpd_crop'));
        end
        
        % We save all images here.
        images{dd} = imread(testFilenameTemp);
    end
    
    %% Plot the camera images.
    PLOTIMAGE = true;
    if (PLOTIMAGE)
        figure;
        figurePosition = [0 0 800 800];
        set(gcf,'position',figurePosition);
        sgtitle(legendHandles{cc});
        
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
        sgtitle(legendHandles{cc});
        
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
            [contrastsRawTemp ESFTemp] = GetImgContrast(images{dd},'minPeakDistance',minPeakDistance);
            contrastsRaw{dd} = contrastsRawTemp;
            meanContrasts{dd} = mean(contrastsRawTemp);
            stdErrorContrasts{dd} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));

            % Print out the Edge spread function to compare across the
            % channels.
            ESF{cc,dd} = ESFTemp;
        end
        
        % Save the plot if you want.
        SAVETHEPLOT = false;
        if (SAVETHEPLOT)
            testFileFormat = '.tiff';
            testFilenameTemp = sprintf('%d nm',peaks_spd(cc));
            saveas(gcf,append(testFilenameTemp,testFileFormat));
            disp('Plot has been saved successfully!');
        end
    end
    
    % Collect the mean contrast results.
    meanContrasts_all(:,cc) = cell2mat(meanContrasts);
end

%% Plot MTF results on one panel.
figure; clf; hold on;
for cc = 1:nTargetChs
    % Printed pattern.
    plot(cell2mat(targetCyclePerDeg),meanContrasts_all(:,cc),...
        '.-','markersize',20);
end

% Get marker color info so that we can use the same one for the same
% channel.
f = flip(get(gca,'children'));

ylim([0 1]);
xlabel('Spatial Frequency (cpd)','fontsize',15);
ylabel('Mean Contrasts','fontsize',15);
xticks(cell2mat(targetCyclePerDeg));
legend(legendHandles,'location','southwest','fontsize',15);

% This is one cycle contrast results. We will calculate the compensated
% contrasts by dividing when plotting the results.
%
% These are the contrats when measuring only one cycle (so, half black on
% the left and the other half as white on the right)
%
% Measured on 0905.
refContrasts = [0.8811    0.8756    0.8938    0.8891    0.8863    0.8870    0.8853    0.8868];

%% Plot MTF comparing with SACCSFA results.
figure; clf;
figureSize = [0 0 1000 500];
set(gcf,'position',figureSize);
sgtitle('MTF comparison: Camera vs. SACCSFA', 'fontsize', 15);
ss = 1;

% Normalize the contrast by dividing the single cycle contrast.
meanContrastsNorm_all = meanContrasts_all./refContrasts;

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
    
    % Set target channel temporarily.
    targetChTemp = idxChComparison(cc);
    % Printed pattern.
    meanContrastsTemp = meanContrastsNorm_all(:,targetChTemp);
    plot(cell2mat(targetCyclePerDeg),meanContrastsTemp,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % SACCSFA.
    if ismember(targetChTemp,idxChComparison)
        plot(cell2mat(targetCyclePerDeg),meanContrasts_SACCSFA(:,ss),...
            'ko-','markeredgecolor','k','markerfacecolor','r','markersize',10);
        ss = ss + 1;
        
        % Add legend.
        legend('Camera','SACCSFA','location','southeast','fontsize',11);
    else
        % Add legend.
        legend('Camera','location','southeast','fontsize',11);
    end
    
    ylim([0 1.1]);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Mean Contrasts','fontsize',15);
    xticks(cell2mat(targetCyclePerDeg));
    title(sprintf('%d nm', peaks_spd_SACCSFA(cc)), 'fontsize', 15);
end

%% Plot the MTF of SACCSFA under assumption using a perfect camera.
%
% Here we divide the MTF by camera MTF.
meanContrastsSACCSFAPerfect = meanContrasts_SACCSFA./meanContrastsNorm_all(:,idxChComparison);

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
for dd = 1:nSFs
    subplot(5,1,dd); hold on;
    
    % Channel.
    for cc = 1:nChComparison
        idxTemp = idxChComparison(cc);
        ESFTemp = ESF{cc,dd};
        plot( (ESFTemp-min(ESFTemp)) ./ (max(ESFTemp)-min(ESFTemp)) );
        
        % Generate texts for the legend.
        legendHandlesESF{cc} = append(num2str(peaks_spd(idxTemp)),' nm');
    end

    % Set each graph in format.
    title(sprintf('%d cpd',targetCyclePerDeg{dd}),'fontsize',15);
    legend(legendHandlesESF,'fontsize',11,'location','southeast','fontsize',10);
    xlabel('Pixel position (horizontal)','fontsize',12);
    ylabel('Normalized dRGB','fontsize',12);
end

%% Doing FFT for further analysis.
%
% Read the signal
numChannel = 2;
SF = 1;
signal = ESF{numChannel,SF};

% Define the parameters of the complex signal
% samplingRate = size(signal,2); 
samplingRate = length(signal);

% Set duration of the signal in seconds. We set arbitrary number here for
% calculation.
duration = 1; 

% Generate the time vector
t = linspace(0,duration,samplingRate);

% Perform the Fourier Transform here.
N = length(signal);
frequencies = (0:N-1) * (samplingRate / N);
signal_fft = fft(signal);
signal_fft_magnitude = abs(signal_fft) / N;

% Find the fundamental frequencies
[~, index1] = max(signal_fft_magnitude);
fundamental_frequency1 = frequencies(index1);

% Find the secound fundamental frequencies. We will remove the peak for the
% first frequency
signal_fft_magnitude(index1) = 0; 
[~, index2] = max(signal_fft_magnitude);
fundamental_frequency2 = frequencies(index2);

% Find the third fundamental frequencies
signal_fft_magnitude(index2) = 0; 
[~, index3] = max(signal_fft_magnitude);
fundamental_frequency3 = frequencies(index3);

% Plot the original signal and its frequency spectrum
figure;
xlim([-1 N]);
subplot(2, 1, 1);
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Target signal');

subplot(2, 1, 2);
plot(frequencies, signal_fft_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum');

fprintf('Fundamental Frequency 1: %.2f Hz\n', fundamental_frequency1);
fprintf('Fundamental Frequency 2: %.2f Hz\n', fundamental_frequency2);
fprintf('Fundamental Frequency 3: %.2f Hz\n', fundamental_frequency3);

% Plot the reconstructrued signal.
% Reconstruct the waveform using the fundamental frequency
reconstructed_signal1 = sin(2 * pi * fundamental_frequency1 * t);
reconstructed_signal2 = sin(2 * pi * fundamental_frequency2 * t);
reconstructed_signal3 = sin(2 * pi * fundamental_frequency3 * t);

% Plot the original and reconstructed signals
figure;
subplot(4, 1, 1);
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');

subplot(4, 1, 2);
plot(t, reconstructed_signal1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal (Fundamental Frequency 1)');

subplot(4, 1, 3);
plot(t, reconstructed_signal2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal (Fundamental Frequency 2)');

subplot(4, 1, 4);
plot(t, reconstructed_signal3);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal (Fundamental Frequency 3)');

%% To use fmincon to fit the signal.
%
% Define the observed waveform data (replace this with your actual data)
x_data = t;
observed_waveform = signal-median(signal);

% Define the objective function for optimization
objective_function = @(params) norm(params(1)*sin(2*pi*params(2)*x_data + params(3)) - observed_waveform);

% Initial guess for parameters a, b, and c
initial_guess = [1, 1/0.2, 0];

% Lower and upper bounds for parameters
lb = [0, 0, -pi];
ub = [Inf, Inf, pi];

% Call fmincon to optimize the parameters
options = optimoptions('fmincon', 'Display', 'iter');
optimized_params = fmincon(objective_function, initial_guess, [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
a_optimized = optimized_params(1);
b_optimized = optimized_params(2);
c_optimized = optimized_params(3);

% Generate the fitted waveform using the optimized parameters
fitted_waveform = a_optimized*sin(2*pi*b_optimized*x_data + c_optimized);

% Plot the observed and fitted waveforms
figure;
plot(x_data, observed_waveform, 'b', 'DisplayName', 'Observed Waveform');
hold on;
plot(x_data, fitted_waveform, 'r', 'DisplayName', 'Fitted Waveform');
xlabel('x');
ylabel('Amplitude');
legend('show');

% Display optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('a = %.4f\n', a_optimized);
fprintf('b = %.4f\n', b_optimized);
fprintf('c = %.4f\n', c_optimized);


% Different way to do fmincon, but lock for now.
%
% Set bounds for parameter x to 0 and 1.
% x0 = [1, 1/(1/15/(2*pi)), 0];
% lb = [0.1, 0.1, -pi];
% ub = [inf, inf, pi];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% options = optimset('fmincon');
% 
% x_found = fmincon(@(x) SineFitSearchFunction(signal, t, x), ...
%     x0, A, b, Aeq, beq, lb, ub, [], options);
% A = x_found(1);
% f = x_found(2);
% phi = x_found(3);
% 
% figure; hold on;
% plot(t,signal);
% plot(t,A*sin(2*pi*f*t+phi));

%% Different method 
figure;
x = linspace(0,1,N)';
y = signal';
plot(x,y,'.')

mdl = fittype('a*sin(b*x+c)+d','indep','x');
fittedmdl2 = fit(x,y,mdl,'start',[rand(),1/(1/5/(2*pi)),rand(),rand()])

figure;
plot(fittedmdl2)
hold on
plot(x,y,'.')
hold off
