% SACC_runExperiment
%
% This is for running a psychphysical experiment for the SACC project. It
% basically contains the three parts - Initialize / Running trials / Close.

% History:
%    01/03/22 smo    Started on it.

%% Initialize.
clear all; close all;

%% Set paratmeters.
%
% If you are not acutally running the experiment, set TURNONSCREEN as false
% and just image will show up instead of displaying on PTB. It will be
% helpful to check and debug.
initialScreenSettings = [1 1 1];
nTrials = 10;

TURNONSCREEN = false;

%% Open the screen.
if(TURNONSCREEN)
    OpenPlainScreen(initialScreenSettings);
end

%% Running trials.
%
% Load the saved data for getting a test image.
% Image is stored in 'ScreenSettingsImage'
%
% Maybe we want to save the test images in the same .mat file so that we
% can load the file and read the images from the same source of the file.
conditionName = 'LminusMSmooth';

% Load test image file.
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theData = load(testFilename);
end
image = theData.screenSettingsImage;

% Display test image.
%
% Now it just displays an image, but this part will be substituted with a
% separate function displaying the image using PTB later on.
for tt = 1:nTrials
    figure; clf;
    imshow(image);
    title(append('Trial ',num2str(tt)),'fontsize',15);
    fprintf('Test image 1 - trial %d is displaying and waiting for the key is pressed... \n',tt);
    
    % Get a response either Yes or No.
    %
    % It can be also done by using 'ginput'. But, here we used the function
    % waitforbuttonpresss.
    %
    % Following is the ASCII allocated number for the keyboards.
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    gettingResponse = waitforbuttonpress;
    response(tt) = double(get(gcf,'CurrentCharacter'));
    close;
    fprintf('     Key input has been received! \n');
end

fprintf('Test image 1 evalaution complete! \n');

%% Close.
if(TURNONSCREEN)
    CloseScreen;
end

% Save the response data.
%
% We may want to change the folder later on. But, now it is set to the same
% folder where the test image is stored.
if(TURNONSCREEN)
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = fullfile(testFiledir,sprintf('expResponseData_%s',conditionName));
        save(testFilename,'response');
    end
end