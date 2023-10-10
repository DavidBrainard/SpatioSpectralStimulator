% SpectralCalCheck_ver3.
%
% This routine reads in the image settings that used in the experiment for
% the SACC proejct, and compare its predicted contrasts with the desired
% contrasts per each pixel of the image.
%
% This is the third version of the validation program, and it does do
% calculation of the settings using the Standard method which was actually
% used in the experiment.
%
% The biggest difference from the previous versions is that this version
% uses the settings from the images that were actually presented in the
% experiment, rather than using the contrast levels to test.
%
% See also:
%    SpectralCalCheck, SpectralCalCheck_ver2

% History:
%    10/05/2023  smo  Started on it.

%% Initialize.
clear; close all;

%% Read in the image settings that used in the experiment.
%
% Set the image type to load either 'normal' or 'high',
imageType = 'normal';
imageSpatialFrequency = 18;

% Load the test image data here.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
    % We set the test file name differently over the image type either
    % 'normal' or 'high'.
    olderDate = 0;
    switch imageType
        case 'normal'
            testFilename = GetMostRecentFileName(testFiledir,sprintf('RunExpData_%d_cpd_',imageSpatialFrequency),'olderDate',olderDate);
        case 'high'
            testFilename = GetMostRecentFileName(testFiledir,'RunExpData_high_18_cpd_','olderDate',olderDate);
    end
    % Load the image date here.
    imageData = load(testFilename);
end

% Get the date of the experiment. We will match this with the validation
% data which will be loaded in the next step.
numExtract = regexp(testFilename,'\d+','match');
yearStr = numExtract{2};
monthStr = numExtract{3};
dayStr = numExtract{4};
dateStrExp = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);

% Load the highest contrast image to test. Test images are saved in an
% ascending order of the contrast, so we will load the last image in the
% array.
testImage = imageData.sceneParamsStruct.predefinedRGBImages{end};

% Display the image if you want.
SHOWTESTIMAGE = false;
if (SHOWTESTIMAGE)
    figure; imshow(testImage);
    title('Following image is going to be tested');
end

% Convert the test image settings to cal format for calculation.
imageTestSettingsCal = ImageToCalFormat(testImage);
nTestSettings = size(imageTestSettingsCal,2);

% Get the image background settings. Here we simply take the very first
% pixel of it which should be the background.
imageBgSettings = imageTestSettingsCal(:,1);

%% Load the image validation data.
%
% Here we load the validation data which contains the measured primaries
% info so that we can calculate the predicted contrasts of the settings.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'CheckCalibration');
    
    % We make a loop here to find the file that matches the date of the
    % experiment. The date of the experiment is in 'dateStrExp' and we will
    % match it with the date of validation which is in 'dateStrVal'.
    %
    % Fix 'olderDate' to 0 here so that it can search from the most recent
    % data.
    olderDate = 0;
    while 1
        % Get test file name.
        testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck_','olderDate',olderDate);
        
        % Extract the date from the file name.
        numExtract = regexp(testFilename,'\d+','match');
        yearStr = numExtract{1};
        monthStr = numExtract{2};
        dayStr = numExtract{3};
        dateStrVal = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);
        
        % Check if two dates match, if it does, stop the loop and load the
        % file.
        if strcmp(dateStrExp,dateStrVal)
            % If we test 'normal' image set, we load the data second to the
            % most recent one. This is because, when we used the 'high'
            % image set, there will be two different validation data with
            % the same date, we always did the validation of 'normal' image
            % set first. Not the most elaborate way, but we do it in this
            % way for now.
            if strcmp(imageType,'normal')
                testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck_','olderDate',olderDate+1);
            end
            fprintf('Target file date: (%s) / Found file date: (%s). Date matched! Data will be loaded. \n',dateStrExp,dateStrVal);
            break;
        else
            % Disable the comment for now. We can re-able it if we want.
            %             fprintf('Target file date: (%s) / Found file date: (%s). We will search again... \n',dateStrExp,dateStrVal);
            olderDate = olderDate+1;
        end
    end
    
    % Load the file here.
    theValData = load(testFilename);
end

% All variables we need are stored in the validation data that we just
% loaded in 'theValData'. Extract some variables for convenience.
%
% Note that 'screenCalObj' is storing the measured primaries in the struct
% 'P_device' as we saved out in that way.
screenCalObj = theValData.screenCalObj;
nPrimaries = size(theValData.targetScreenSpdMeasured,2);

%% Calculate the predicted contrasts of the test image presented in the experiment.
%
% Calculate the cone excitations of the backgrond settings for calculation
% of the contrasts.
imageBgPrimaries = SettingsToPrimary(screenCalObj,imageBgSettings);
imageBgExcitations = PrimaryToSensor(screenCalObj,imageBgPrimaries);

% Get contrasts of the test settings per each pixel.
imageTestPrimaries = SettingsToPrimary(screenCalObj,imageTestSettingsCal);
imageTestExcitations = PrimaryToSensor(screenCalObj,imageTestPrimaries);
imageTestContrastsCal = ExcitationsToContrast(imageTestExcitations,imageBgExcitations);

%% Calculate the desired contrasts of the test image per each pixel.
%
% Here, we will calculate the desired contrast of each pixel on the test
% image so that we can compare it with the actual contrast pixel by pixel.
%
% We use the same sub-functions here that were used to make the test images
% in the experiment.
%
% Get screen size object. Same as the test image.
[~,screenSizeObject,~] = SetupISETBioDisplayObject(imageData.colorDirectionParams,screenCalObj,'verbose',false);

% Read out some variables from the image data to make monochrome gabor. We
% will read out from the image data.
stimulusSizeDeg = imageData.spatialTemporalParams.stimulusSizeDeg;
gaborSdDeg = imageData.spatialTemporalParams.gaborSdDeg;
sineFreqCyclesPerDeg = imageData.spatialTemporalParams.sineFreqCyclesPerDeg;
sineImagePhaseShiftDeg = imageData.spatialTemporalParams.sineImagePhaseShiftDeg;

% Make monochrome gabor image here.
nQuantizeBits = 14;
[MonochromeGaborImage, ~, rawMonochromeContrastGaborCal,~,~,~,~] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,...
    'sineImagePhaseShiftDeg',sineImagePhaseShiftDeg,'verbose',false,'nQuantizeBits',nQuantizeBits);

% Plot the monochrome image if you want.
SHOWMONOGABORIMAGE = false;
if (SHOWMONOGABORIMAGE)
    figure;
    imshow(MonochromeGaborImage{1});
    title('Monochrome gabor image');
end

% Calculate the desired contrast for all pixels using the monochrome gabor
% image we made from the above. We muliply the target contrast and target
% stimulus contrast direction (L-M) to the monochrome gabor image.
desiredContrastGaborCal = imageData.colorDirectionParams.spatialGaborTargetContrast * imageData.colorDirectionParams.targetStimulusContrastDir * cell2mat(rawMonochromeContrastGaborCal);

% Calculate the image contrast. We will set the L-M contrast as a positive
% sign and -(L-M) as a negative sign.
for cc = 1:size(desiredContrastGaborCal,2)
    desiredImageContrastGaborCalTemp = sqrt(sum(desiredContrastGaborCal(:,cc).^2));
    if desiredContrastGaborCal(1,cc) > desiredContrastGaborCal(2,cc)
        % L-M
        desiredImageContrastGaborCal(cc) = desiredImageContrastGaborCalTemp;
    else
        % -(L-M)
        desiredImageContrastGaborCal(cc) = -desiredImageContrastGaborCalTemp;
    end
end

%% Plot the results: Predicted cone contrast vs. Desired cone contrast.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};
sgtitle(sprintf('Desired cone contrasts vs. Predicted cone contrasts in the image (N=%d)',nTestSettings));

% Make a loop for plotting each cone case.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % We will plot Standard and PointCloud methods separately.
    plot(desiredContrastGaborCal(pp,:),imageTestContrastsCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    plot(desiredContrastGaborCal(pp,:),desiredContrastGaborCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired cone contrast','fontsize',15);
    ylabel('Predicted cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('Test Image','Desired','location','southeast','fontsize',13);
end

%% Plot the results: Predicted cone contrast vs. Desired image contrast.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};
sgtitle(sprintf('Desired image contrast vs. Actual cone contrasts in the image (N=%d)',nTestSettings));

% Sort the image contrast in an ascending order.
[desiredImageContrastGaborCalSorted I] = sort(desiredImageContrastGaborCal,'ascend');

% Sort the cone contrasts in the same order as well.
imageTestContrastsCalSorted = imageTestContrastsCal(:,I);
desiredContrastGaborCalSorted = desiredContrastGaborCal(:,I);

% Marginal contrast found from the validation of the nominal contrasts.
% This was found from the routine SpectralCalCheck_ver2.
%
% Normal image set. Put negative sign to the image contrast, because we set
% (L-M) as a positive sign and -(L-M) as a negative sign.
coneContrastNormal = [-0.0427 0.0369 -0.0028];
imageContrastNormal = -sqrt(sum(coneContrastNormal.^2));

% High image set.
coneContrastHigh = [-0.0579 0.0590 -0.0003];
imageContrastHigh = -sqrt(sum(coneContrastHigh.^2));

% Make a loop for plotting each cone case.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % We will plot Standard and PointCloud methods separately.
    plot(desiredImageContrastGaborCalSorted,imageTestContrastsCalSorted(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    plot(desiredImageContrastGaborCalSorted,desiredContrastGaborCalSorted(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    
    % Mark the cut-off contrast point if you want.
    TESTONEPOINT = true;
    if(TESTONEPOINT)
        % Set the index within the test contrast array to plot.
        %
        % For 12 cpd, index = 290.
        switch imageType
            case 'normal'
                index = 1850;
            case 'high'
                index = 1500;
        end
        
        % Contrast cut-off point.
        plot(desiredImageContrastGaborCalSorted(index),imageTestContrastsCalSorted(pp,index),'o','MarkerSize',14,'markerfacecolor','y','markeredgecolor','k');
        % Contrast cut-off line.
        plot(ones(1,2)*desiredImageContrastGaborCalSorted(index),[-0.1 0.1],'color',[1 1 0 0.7],'linewidth',5);
        
        % Calculate the marginal contrast from the nominal contrasts to
        % print out.
        marginalContrast = sqrt(sum(imageTestContrastsCalSorted(:,index).^2));
        % We will only print once by setting 'pp' to 1.
        if pp == 1
            fprintf('Good maximum image contrast = (%.4f / %.4f) / Log sensitivity = (%.4f) \n',marginalContrast, abs(desiredImageContrastGaborCalSorted(index)), log10(1/marginalContrast));
        end
    end
    
    % Add the marginal contrasts found from the earlier testing.
    switch imageType
        case 'normal'
            plot(imageContrastNormal,coneContrastNormal(pp),'o','markersize',10,'markerfacecolor','k','markeredgecolor','k');
        case 'high'
            plot(imageContrastHigh,coneContrastHigh(pp),'o','markersize',10,'markerfacecolor','k','markeredgecolor','k');
    end
    
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired image contrast','fontsize',15);
    ylabel('Predicted cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    grid on;
    legend('Test Image','Desired','location','southeast','fontsize',13);
end
