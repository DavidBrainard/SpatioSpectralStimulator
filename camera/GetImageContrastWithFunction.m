% GetImageContrastWithFunction.
%
% It calculates the image contrasts using the function.

% History:
%    06/13/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Load images.
%
% Set target spatial frequency.
targetCyclePerDeg = {3, 6, 9, 12, 18};

% Load image here.
nImages = length(targetCyclePerDeg);
for tt = 1:nImages
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
        testFiledir = fullfile(testFiledir,'Camera','ChromaticAberration');
        testFilename = GetMostRecentFileName(testFiledir,append(num2str(targetCyclePerDeg{tt}),'cpd_crop'));
        image{tt} = imread(testFilename);
    else
        error('Cannot find data file');
    end
end

% Calculate contrast.
for tt = 1:nImages
    subplot(5,1,tt);
    contrast(tt) = GetImgContrast(image{tt});
end