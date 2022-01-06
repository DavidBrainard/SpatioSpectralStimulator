% SACC_runExperiment
%
% This is for running a psychphysical experiment for the SACC project. It
% basically contains the three parts - Initialize / Running trials / Close.

% History:
%    01/03/22 smo    Started on it.

%% Initialize.
clear; close all;

%% Set paratmeters.
%
% If you are not acutally running the experiment, set TURNONSCREEN as false
% and just image will show up instead of displaying on PTB. It will be
% helpful to check and debug.
initialScreenSettings = [1 1 1];
nTestImages = 2;
nTrials = 3;

VERBOSE = true;
TURNONSCREEN = false;

%% Open the screen.
if (TURNONSCREEN)
    OpenPlainScreen(initialScreenSettings,'verbose',VERBOSE);
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
for ii = 1:nTestImages
    for tt = 1:nTrials
        figure; clf;
        imshow(image);
        title(append('Test Image ',num2str(ii),' - Trial ',num2str(tt)),'fontsize',15);
        if (VERBOSE)
            fprintf('Test image %d - trial %d is displaying and waiting for the key is pressed... \n',ii,tt);
        end
        
        % Get a response either Yes or No.
        %
        % It can be also done by using 'ginput'. But, here we used the function
        % waitforbuttonpresss.
        %
        % Following is the ASCII allocated number for the keyboards.
        %
        % 28 leftarrow
        % 29 rightarrow
        % 30 uparrow
        % 31 downarrow
        leftArrow  = 28;
        rightArrow = 29;
        gettingResponse = waitforbuttonpress;
        response(tt,ii) = double(get(gcf,'CurrentCharacter'));
        close;
        
        if (VERBOSE)
            fprintf('     Key input has been received! \n');
        end
    end
    if(VERBOSE)
        fprintf('Test image %d evalaution complete! \n',ii);
    end
end

%% Show the results.
%
% Convert the response into 0 / 1
response(response == leftArrow)  = 0;
response(response == rightArrow) = 1;

for ii = 1:nTestImages
    propCorrect(:,ii) = sum(response(:,ii)==1)/nTrials;
end

% Plot it.
figure; clf;
plot(linspace(1,nTestImages,nTestImages),propCorrect,'r.','markersize',15);
xlabel('Test image','fontsize',15);
ylabel('Proportion correct (%)','fontsize',15);
xlim([0 nTestImages+1]);
ylim([0 1]);

% Show only integer value on x axis
curTick = get(gca, 'xTick');
xticks(unique(round(curTick)));

%% Close.
if (TURNONSCREEN)
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
end