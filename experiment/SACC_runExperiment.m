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
timeDelayBtwImages = 1;

VERBOSE = true;
TURNONSCREEN = true;

%% Load the test images for the experiment.
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

%% Running trials using PTB.
%
% This part will be used for the actual experiment displaying the test
% image on the projector using PTB.
if (TURNONSCREEN)
    % Open the screen ready.
    [window windowRect] = OpenPlainScreen(initialScreenSettings,'verbose',VERBOSE);
    
    % Display the test images here.
    for ii = 1:nTestImages
        for tt = 1:nTrials
            % First Image.
            SetScreenImage(image, window, windowRect,'verbose',VERBOSE);
%             MakeBeepSound;
            
            % Make a time delay before displaying the other image of the
            % pair.
            for dd = 1:timeDelayBtwImages;
                pause(1);
            end
            
            % Second Image.
            SetScreenImage(image, window, windowRect,'verbose',VERBOSE);
%             MakeBeepSound;
            
            if (VERBOSE)
                fprintf('Test image %d - trial %d is displaying and waiting for the key is pressed... \n',ii,tt);
            end
            
            % Get a key stroke response here.
            gettingResponse = waitforbuttonpress;
            response(tt,ii) = double(get(gcf,'CurrentCharacter'));
            close all;
            if (VERBOSE)
                fprintf('     Key input has been received! \n');
            end
        end
        if(VERBOSE)
            fprintf('Test image %d evalaution complete! \n',ii);
        end
    end
end

%% Running trials - Simulation without using PTB.
%
% This part does not use PTB and just diplay test images side-by-side on
% Figure and record key stroke responses, which would be helpful to check
% the sequence of the experiment before running it with the patients.
if (~TURNONSCREEN)
    % Resize the image as desired here.
    imageMagnificationFactor = 1.5;
    image = imresize(image,imageMagnificationFactor);
    imageSize = size(image);
    imageXPixel = imageSize(1);
    imageYPixel = imageSize(2);
    
    % Get the display info here.
    screenSize = get(0,'screensize');
    screenXPixel = screenSize(3);
    screenYPixel = screenSize(4);
    
    % Display test image.
    %
    % Now it just displays an image, but this part will be substituted with a
    % separate function displaying the image using PTB later on.
    for ii = 1:nTestImages
        for tt = 1:nTrials
            figure; clf;
            
            % Set the position and the size of the test image.
            imageFig = figure;
            x = screenXPixel*0.2;
            y = screenYPixel*0.2;
            width = imageXPixel*2;
            height = imageYPixel;
            set(imageFig, 'Position', [x y width height])
            
            % Left side image.
            subplot(1,2,1); imshow(image);
            title(append('Test Image ',num2str(ii),' - Trial ',num2str(tt)),'fontsize',15);
            
            % Right side image.
            subplot(1,2,2); imshow(image);
            
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
            gettingResponse = waitforbuttonpress;
            response(tt,ii) = double(get(gcf,'CurrentCharacter'));
            close all;
            
            if (VERBOSE)
                fprintf('     Key input has been received! \n');
            end
        end
        if(VERBOSE)
            fprintf('Test image %d evalaution complete! \n',ii);
        end
    end
end

%% Show the results.
%
% Convert the response into 0 / 1
leftArrow  = 28;
rightArrow = 29;
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