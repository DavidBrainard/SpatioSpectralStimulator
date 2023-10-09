% DisplayGaborPatch. 
%
% This displays an image on the DLP using the Psychtoolbox. 
%
% History:
%    10/09/23  smo    - Modified it.

%% Initialize.
close all; clear;

%% Load the image data.
imageType = 'normal';
olderDate = 5;
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
    % We set the test file name differently over the image type either
    % normla or high.
    switch imageType
        case 'normal'
            testFilename = GetMostRecentFileName(testFiledir,'RunExpData_18_cpd_','olderDate',olderDate);
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

% Set primary channel settings as we did in the experiment.
channelSettings = imageData.experimentParams.screenPrimarySettings;
SetChannelSettings(channelSettings);

%% PTB pre-setup.
PsychDefaultSetup(2); 

screens = Screen('Screens');    
screenNumber = max(screens);

white = WhiteIndex(screenNumber); 
black = BlackIndex(screenNumber); 

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Read the resolution of the display

%% Display a background as plain screen
%
% Backgroundcolor can be set, but here it is set as the same size between
% the background and the display window, so it will not be shown on the final result
backgroundColor = [255 255 255]; % Range 0-255 (8 bits)
[window, backgroundRect] = Screen('OpenWindow', screenNumber, backgroundColor, windowRect); % Size of the background = Size of the window

%% Display the gabor patch as 'texture'
gaborPatchImage_texture = Screen('MakeTexture', window, testImage); % Set the gaborPatchImage as 'texture' in PTB
Screen('DrawTexture', window, gaborPatchImage_texture, [], backgroundRect);
Screen('Flip', window);
