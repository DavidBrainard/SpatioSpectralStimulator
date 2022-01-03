% SACC_runExperiment
% 
% This is for running a psychphysical experiment for the SACC project. It
% basically contains the three parts - Initialize / Running trials / Close.

% History:
%    01/03/22 smo    Started on it.

%% Initialize.
clear all; close all;

%% Set paratmeters.


%% Open the screen.
initialScreenSettings = [1 1 1];
OpenPlainScreen(initialScreenSettings);

%% Running trials.
%
% Load the saved data for getting a test image.
% Image is stored in 'ScreenSettingsImage'
conditionName = 'LminusMSmooth';

if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theData = load(testFilename);
end

% Display the image.
%
% Now it just displays an image, but this part will be substituted with a
% separate function displaying the image using PTB later on.
image = theData.screenSettingsImage;
imshow(image);

% Get a response either Yes or No.


% Save the response.


%% Close.
CloseScreen;
