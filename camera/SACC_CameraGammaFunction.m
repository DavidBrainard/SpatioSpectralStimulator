% SACC_CameraGammaFunction.
%
% This code is to calculate camera gamma function.
%
% See also:
%    SACC_GetCameraImageOverExposure.

% History:
%    08/14/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
measureDate = '0814';
testChannel = 'Ch7';

%% Load all images here.
digitalInputOptions = [0.1:0.1:1.0];
nImages = length(digitalInputOptions);
for ee = 1:nImages
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Camera','CameraGammaFunction',testChannel);
        digitalInputTemp = digitalInputOptions(ee);
        testFilename = GetMostRecentFileName(testFiledir, sprintf('%s_%.1f%s',testChannel,digitalInputTemp,'crop'));
    else
        error('Cannot find data file');
    end
    % We save all images here. The array looks like {dataType,
    % channel, SF}.
    images{ee} = imread(testFilename);
end

%% Plot the camera images.
PLOTIMAGE = true;
if (PLOTIMAGE)
    figure;
    figurePosition = [0 0 800 800];
    set(gcf,'position',figurePosition);
    
    for ee = 1:nImages
        subplot(5,2, ee);
        imshow(images{ee});
        digitalInputTemp = digitalInputOptions(ee);
        title(sprintf('Image digital input = %.1f',digitalInputTemp),'FontSize',15);
    end
end

%% Plot the mean power over exposure level.
%
% Get mean power of the images.
for ee = 1:nImages
    meanPowerImages(ee) = mean(images{ee},'all');
end
figure; hold on;
xlabel('Digital input (0-1)','fontsize',15);
ylabel('Mean power (dRGB)','fontsize',15);
plot(digitalInputOptions,meanPowerImages,'ro','markersize',11,'markerfacecolor','r','markeredgecolor','k');
