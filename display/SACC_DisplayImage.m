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
    spatialFrequencyOptions = [3 6 9 12 18];
    while 1
        spatialFrequency = input('Which spatial frequency to display? [3, 6, 9, 12, 18] ');
        if ismember(spatialFrequency, spatialFrequencyOptions)
            break
        end
        disp('Choose one of the above spatial frequencies!');
    end
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

% Rotate the image if you want. Here we can rotate the image if you want.
% For SACC experiment, we did use 45-deg rotated image either up left or up
% right.
ROTATEIMAGE = true;
if (ROTATEIMAGE)
    % Get some info on background settings.
    settingsBG = squeeze(testImage(1,1,:));
    settingsBlk = 0;
    
    % Set the rotation direction here either to left or right. 
    imgRotatationDir = 'left'; 
    switch imgRotatationDir
        case 'left'
            imgRotationDeg = 45;
        case 'right'
            imgRotationDeg = -45;
    end
    
    % Image rotation happens here.
    testImage = imrotate(testImage,imgRotationDeg,'crop');
    nPrimaries = size(testImage,3);
    for pp = 1:nPrimaries
        testImage(:,:,pp) = changem(testImage(:,:,pp), settingsBG(pp), settingsBlk);
    end
end

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
