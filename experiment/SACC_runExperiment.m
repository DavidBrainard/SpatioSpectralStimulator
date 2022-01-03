% SACC_runExperiment
%
% This is for running a psychphysical experiment for the SACC project. It
% basically contains the three parts - Initialize / Running trials / Close.

% History:
%    01/03/22 smo    Started on it.

%% Initialize.
clear all; close all;

%% Set paratmeters.
initialScreenSettings = [1 1 1];
nTestImage = 20;

%% Open the screen.
OpenPlainScreen(initialScreenSettings);

%% Running trials.
%
% Load the saved data for getting a test image.
% Image is stored in 'ScreenSettingsImage'
%
% Maybe we want to save the test images in the same .mat file so that we
% can load the file and read the images from the same source of the file. 
conditionName = 'LminusMSmooth';

if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theData = load(testFilename);
end
image = theData.screenSettingsImage;

% Display the test image.
%
% Now it just displays an image, but this part will be substituted with a
% separate function displaying the image using PTB later on.
figure; clf;
imshow(image);
disp('Test image is displaying!');

% Get a response either Yes or No.
disp('Waiting for the key is pressed');

% One way to do it
%
% 28 leftarrow
% 29 rightarrow
% 30 uparrow
% 31 downarrow
gettingResponse = waitforbuttonpress;
response = double(get(gcf,'CurrentCharacter'));
close;

% Another way to do it.
[~,~,response2] = ginput(1);
switch response2
    case 28
    case 29
end
close;

%% Close.
CloseScreen;

% Save the response data.
% 
% We may want to change the folder later on. But, now it is set to the same
% folder where the test image is stored.
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('expResponseData_%s',conditionName));
    save(testFilename,'response');
end