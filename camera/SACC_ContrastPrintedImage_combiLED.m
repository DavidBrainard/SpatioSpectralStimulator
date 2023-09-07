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
            contrastsRawTemp = GetImgContrast(images{dd},'minPeakDistance',minPeakDistance);
            contrastsRaw{dd} = contrastsRawTemp;
            meanContrasts{dd} = mean(contrastsRawTemp);
            stdErrorContrasts{dd} = std(contrastsRawTemp)/sqrt(length(contrastsRawTemp));
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
figureSize = [0 0 1200 500];
set(gcf,'position',figureSize);
sgtitle('MTF comparison: Camera vs. SACCSFA', 'fontsize', 15);
ss = 1;

for cc = 1:nTargetChs
    subplot(2,4,cc); hold on;
    
    % Printed pattern.
    %
    % You can plot either raw contrast or normalized contrast.
    normContrast = true;
    if (normContrast)
        meanContrastsTemp = meanContrasts_all(:,cc)/refContrasts(cc);
    else
        meanContrastsTemp = meanContrasts_all(:,cc);
    end
    plot(cell2mat(targetCyclePerDeg),meanContrastsTemp,...
        'ko-','markeredgecolor','k','markerfacecolor','b', 'markersize',10);
    
    % SACCSFA. Its peak wavelengths were 422, 476, 530, 592, 658 nm.
    idxChComparison = [2 3 5 6 8];
    if ismember(cc,idxChComparison)
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
    title(sprintf('%d nm', peaks_spd(cc)), 'fontsize', 15);
end
