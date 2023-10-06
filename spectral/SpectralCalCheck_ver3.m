% SpectralCalCheck_ver3.
%
% Read in the image settings that used in the experiment for the SACC
% proejct, and check if the image was generated fine.
%
% This is the third version of if, and it does do calculation of the
% settings using Standard method which was actually used in the experiment.
%
% The biggest difference from the previous versions is that this version
% uses the settings from the images that were actually presented in the
% experiment, rather than setting contrast levels to test.
%
% See also:
%    SpectralCalCheck, SpectralCalCheck_ver2

% History:
%    10/05/2023  smo  Started on it.

%% Initialize.
clear; close all;

%% Read in the image settings.
olderDate = 5;
imageType = 'normal';
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
    switch imageType
        case 'normal'
            testFilename = GetMostRecentFileName(testFiledir,'RunExpData_18_cpd_','olderDate',olderDate);
        case 'high'
            testFilename = GetMostRecentFileName(testFiledir,'RunExpData_high_18_cpd_','olderDate',olderDate);
    end
    imageData = load(testFilename);
end

% Get the date of the experiment. We will match this with the validation
numExtract = regexp(testFilename,'\d+','match');
yearStr = numExtract{2};
monthStr = numExtract{3};
dayStr = numExtract{4};
dateStrImage = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);

% Get the highest test image.
testImage = imageData.sceneParamsStruct.predefinedRGBImages{end};

% Display the image if you want.
SHOWTESTIMAGE = true;
if (SHOWTESTIMAGE)
    figure; imshow(testImage);
    title('Following image is going to be tested');
end

% Convert the test image to cal format for calculation.
imageTestSettingsCal = ImageToCalFormat(testImage);

% Get the image background settings.
imageBgSettings = imageTestSettingsCal(:,1);

%% Load the calibration (validation) data.
%
% Set which data to load. Set 'olderDate' to 0 will load the most recent
% data.
olderDate = 0;

% Load the image file.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'CheckCalibration');
    
    % We make a loop here to load the file that matches the above image
    % settings. The validation file has the same date as the image was
    % made.
    while 1
        % Get test file name.
        testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck_','olderDate',olderDate);
        
        % Extract the date from the file name.
        numExtract = regexp(testFilename,'\d+','match');
        yearStr = numExtract{1};
        monthStr = numExtract{2};
        dayStr = numExtract{3};
        dateStrVal = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);
        
        % Check if two dates match, if it does, stop the loop.
        if strcmp(dateStrImage,dateStrVal)
            disp('Image file and validation file dates matched! Validation data will be loaded');
            break;
        else
            disp('Image and Validation file date not matched, we will search it again...');
            olderDate = olderDate+1;
        end
    end
    
    % Load the file.
    theValData = load(testFilename);
end
theData = theValData.theData;

%% Set variables here.
%
% All variables we need are stored in the validation data we just loaded in
% 'theValData'. Extract some variables for convenience.
%
% Load the measured primaries. In the validation data, 'theValData', we
% saved the screenCalObj by substituting its primaries with measured
% primaries, so we will load it from there.
screenCalObj = theValData.screenCalObj;
targetScreenSpdMeasured = screenCalObj.get('P_device');

% Some more variables.
S = theData.S;
wls = SToWls(S);
nPrimaries = size(theData.screenPrimarySpd,2);
nChannels = theData.channelCalObjs{1}.get('nDevices');
channelNInputLevels = size(theData.channelCalObjs{1}.get('gammaInput'),1);
T_cones = theData.T_cones;
spatialGaborTargetContrast = theData.spatialGaborTargetContrast;
targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
targetLambda = theData.targetLambda;

% Control the messages and plots.
verbose = true;

%% Calculate the actual contrasts presented in the test image.
%
% Get cone excitations of the backgrond settings for calculation of the
% contrasts.
imageBgPrimaries = SettingsToPrimary(screenCalObj,imageBgSettings);
imageBgExcitations = PrimaryToSensor(screenCalObj,imageBgPrimaries);

% Get contrasts of the test settings.
imageTestPrimaries = SettingsToPrimary(screenCalObj,imageTestSettingsCal);
imageTestExcitations = PrimaryToSensor(screenCalObj,imageTestPrimaries);
imageTestContrastsCal = ExcitationsToContrast(imageTestExcitations,imageBgExcitations);

%% Calculate the desired contrasts of the test image.
%
% Here, we will calculate the desired contrast of each pixel on the test
% image so that we can compare it with the actual contrast pixel by pixel.
%
% We use the same sub-functions here that were used to make the test images
% in the experiment.
%
% Get screen size object.
[~,screenSizeObject,~] = SetupISETBioDisplayObject(imageData.colorDirectionParams,screenCalObj,'verbose',false);

% Make monochrome gabor image.
stimulusSizeDeg = imageData.spatialTemporalParams.stimulusSizeDeg;
gaborSdDeg = imageData.spatialTemporalParams.gaborSdDeg;
sineFreqCyclesPerDeg = imageData.spatialTemporalParams.sineFreqCyclesPerDeg;
sineImagePhaseShiftDeg = imageData.spatialTemporalParams.sineImagePhaseShiftDeg;
nQuantizeBits = 14;

[rawMonochromeUnquantizedContrastGaborImage, ~, rawMonochromeContrastGaborCal,~,~,~,~] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObject,...
    'sineImagePhaseShiftDeg',sineImagePhaseShiftDeg,'verbose',false,'nQuantizeBits',nQuantizeBits);

% Plot the monochrome image if you want.
SHOWMONOGABORIMAGE = false;
if (SHOWMONOGABORIMAGE)
    figure;
    imshow(rawMonochromeUnquantizedContrastGaborImage{1});
    title('Monochrome gabor image');
end

% Calculate the desired contrast for all pixels.
desiredContrastGaborCal = imageData.colorDirectionParams.spatialGaborTargetContrast * imageData.colorDirectionParams.targetStimulusContrastDir * cell2mat(rawMonochromeContrastGaborCal);
desiredImageContrastGaborCal = sqrt(sum(desiredContrastGaborCal.^2));

%% Compare Desired vs. Actual cone contrasts.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% Make a loop for plotting each cone case.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % We will plot Standard and PointCloud methods separately.
    plot(desiredContrastGaborCal(pp,:),imageTestContrastsCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    plot(desiredContrastGaborCal(pp,:),desiredContrastGaborCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired cone contrast','fontsize',15);
    ylabel('Nominal cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('Test Image','Desired','location','southeast','fontsize',13);
end
