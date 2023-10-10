% SACC_DisplayImage.
%
% This displays an image on the DLP using the Psychtoolbox.
%
% History:
%    10/09/23  smo    - Modified it.

%% Initialize.
close all; clear;

%% Load the image data.
%
% Choose the image type either 'normal' or 'high'.
imageType = 'normal';

% Select spatial frequency to display if we display 'normal' image. For
% high image, it will be fixed to 18 cpd.
if strcmp(imageType,'normal')
    % Choose SF to display here.
    spatialFrequency = 18;
else
    spatialFrequency = 18;
end
olderDate = 5;
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
    % We set the test file name differently over the image type either
    % normla or high.
    switch imageType
        case 'normal'
            testFilename = GetMostRecentFileName(testFiledir,sprintf('RunExpData_%d_cpd_',spatialFrequency),'olderDate',olderDate);
        case 'high'
            testFilename = GetMostRecentFileName(testFiledir,'RunExpData_high_18_cpd_','olderDate',olderDate);
    end
    % Load the image date here.
    imageData = load(testFilename);
end

% Load the highest test image to test. Test images are saved in an
% ascending order of the contrast, so we will load the last image in the
% array.
testImage = imageData.sceneParamsStruct.predefinedRGBImages{end};

% Load the primary channel settings we used in the experiment.
channelSettings = imageData.experimentParams.screenPrimarySettings;

%% Open the projector and set the channel settings.
%
% Turn on the projector and display anything to start.
initialScreenSetting = [0 0 0]';
[window windowRect] = OpenPlainScreen(initialScreenSetting);

% Set channel settings here. We loaded the settings from the saved data.
SetChannelSettings(channelSettings);

%% Make PTB image texture.
%
% We will use the same function that we used in the experiment to display
% the image.
[imageTexture imageWindowRect rng] = MakeImageTexture(testImage, window, windowRect, ...
    'addNoiseToImage', true, 'verbose', false);

%% Flip the PTB texture to display the image on the projector.
FlipImageTexture(imageTexture, window, imageWindowRect,'verbose',false);
fprintf('Image is now displaying... Current image = (%d cpd) \n',spatialFrequency);
