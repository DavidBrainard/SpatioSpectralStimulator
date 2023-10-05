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
olderDate = 0;
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
    testFilename = GetMostRecentFileName(testFiledir,'RunExpData_18_cpd','olderDate',olderDate);
    imageData = load(testFilename);
end

% Get the date of the experiment. We will match this with the validation
numExtract = regexp(testFileNameContrast,'\d+','match');
yearStr = numExtract{3};
monthStr = numExtract{4};
dayStr = numExtract{5};
dateStr = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);


% Get the highest test image.
testImage = imageData.sceneParamsStruct.predefinedRGBImages{end};

% Display the image if you want.
figure; imshow(testImage);
title('Following image is going to be tested');

% Extract the unique image settings from the above image.
for xx = 1:size(testImage,1)
    for yy = 1:size(testImage,2)
        for cc = 1:size(testImage,3)
            testImageSettingsCal(cc,yy+(xx-1)*size(testImage,2)) = testImage(xx,yy,cc);
        end
    end
end

% Get background settings.
bgSettings = testImageSettingsCal(:,1);

% Find the unique test settings. We will find any settings different from
% background.
nSettings = size(testImageSettingsCal,2);
testSettings = [];
for tt = 1:nSettings
    testSettingTemp = testImageSettingsCal(:,tt);
    % Sort out the background settings.
    if (testSettingTemp(1) == bgSettings(1) && testSettingTemp(2) == bgSettings(2) && testSettingTemp(3) == bgSettings(3))
        continue;
    else
        % Here, we collect the test settings having different settings from
        % the background. 'testSettings' contain all settings presented on
        % the test image in all pixels.
        testSettings(:,end+1) = testSettingTemp;
    end
end

%% Load the calibration (validation) data.
olderDate = 0;

% Load the image file.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'CheckCalibration');
    testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck_','olderDate',olderDate+1);
    theData = load(testFilename);
end

% Match the data format to run it smoothly.
screenCalObj_temp = theData.screenCalObj;
theData = theData.theData;
theData.screenCalObj = screenCalObj_temp;

%% Set variables here.
%
% All variables we need are stored in the image data we just loaded in
% 'theData'. Extract some for convenience.
targetScreenSpd = theData.screenCalObj.get('P_device');
S = theData.S;
wls = SToWls(S);
nPrimaries = size(theData.screenPrimarySpd,2);
nChannels = theData.channelCalObjs{1}.get('nDevices');
channelNInputLevels = size(theData.channelCalObjs{1}.get('gammaInput'),1);
T_cones = theData.T_cones;

% Control the messages and plots.
verbose = true;

if isfield(theData,{'spatialGaborTargetContrast','targetScreenPrimaryContrast','targetLambda'})
    spatialGaborTargetContrast = theData.spatialGaborTargetContrast;
    targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
    targetLambda = theData.targetLambda;
end

% Load the measured primaries.
targetScreenSpdMeasured = theData.screenCalObj.cal.P_device;

%% Get target settings from the test image..
%
% THIS PART WILL BE SUBSTITUTED WITH THE LOADED IMAGE SETTINGS FROM THE EXPERIMENT.
% Here, we calculate background excitation and excitation of each pixel.
%
% Calculate the desired cone cone contrasts.
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
desiredImageContrastCheckCal = theData.spatialGaborTargetContrast * rawMonochromeContrastCheckCal;
desiredContrastCheckCal = theData.spatialGaborTargetContrast * theData.targetStimulusContrastDir * rawMonochromeContrastCheckCal;
desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,screenBgExcitations);

%% Here we calculate the settings to achieve the desired contrast.
%
% Get primaries using standard calibration code, and desired spd without
% quantizing.
standardPrimariesGaborCal = SensorToPrimary(screenCalObj,desiredExcitationsCheckCal);
standardDesiredSpdGaborCal = PrimaryToSpd(screenCalObj,standardPrimariesGaborCal);

% Gamma correct and quantize. Then convert back from the gamma
% corrected settings.
standardSettingsGaborCal = PrimaryToSettings(screenCalObj,standardPrimariesGaborCal);
standardPredictedPrimariesGaborCal = SettingsToPrimary(screenCalObj,standardSettingsGaborCal);
standardPredictedExcitationsGaborCal = PrimaryToSensor(screenCalObj,standardPredictedPrimariesGaborCal);
standardPredictedContrastGaborCal = ExcitationsToContrast(standardPredictedExcitationsGaborCal,screenBgExcitations);

%% Compare Desired vs. Actual cone contrasts.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
sgtitle('Contrast: Standard vs. Point cloud methods');
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% Make a loop for plotting each cone case. We decided to plot the results
% separately over the calculation methods, Standard and PointCloud. Making
% a loop in this way is not the most elaborate, so we may want to update
% this part later on.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % Save the calculation type per each order.
    calculationMethod = 'Standard';
    
    % We will plot Standard and PointCloud methods separately.
    plot(desiredContrastCheckCal(pp,:),standardPredictedContrastGaborCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired cone contrast','fontsize',15);
    ylabel('Nominal cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('Standard','location','southeast','fontsize',13);
end

%% Compare actual cone contrast over desired test image contrast.
% After the meeting (Semin and David) on 09/26/23, we added this part.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
sgtitle('Desired image contrast vs. Nominal contrast');
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% Again, not fancy loop. Same graph as the above, but different way to look
% at it with updated x-axis as the desired image contrast.
nPlots = nPrimaries*2;
for pp = 1:nPlots
    subplot(2,3,pp); hold on;
    
    % Save the calculation type per each order.
    if ismember(pp,[1 2 3])
        calculationMethod = 'Standard';
    elseif ismember (pp, [4 5 6])
        calculationMethod = 'PointCloud';
    end
    
    % We will plot Standard and PointCloud methods separately.
    if pp <= nPrimaries
        % Standard method.
        plot(desiredImageContrastCheckCal,standardPredictedContrastGaborCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
        plot(desiredImageContrastCheckCal,desiredContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
        
        % Search for a marginal contrast making a good image. Set the
        % index from 1 to 202 (number of nominal contrasts to test) so that
        % you can visually check from when the contrast reproduction goes
        % bad.
        %
        % For now, we set the cone contrasts within the same index of the
        % test contrats.
        %
        % We will mark the cut-off contrast points as a line and a point on
        % the figure.
        TESTONEPOINT = true;
        if(TESTONEPOINT)
            % We found index = 183 for Normal image set, and index = 185
            % for high image set.
            index = 183;
            % Contrast cut-off point.
            plot(desiredImageContrastCheckCal(index),standardPredictedContrastGaborCal(pp,index),'o','MarkerSize',14,'markerfacecolor','y','markeredgecolor','k');
            % Contrast cut-off line.
            plot(ones(1,2)*desiredImageContrastCheckCal(index),[-0.1 0.1],'color',[1 1 0 0.7],'linewidth',5);
            
            % Calculate the marginal contrast from the nominal contrasts to
            % print out.
            marginalContrast = sqrt(sum(standardPredictedContrastGaborCal(:,index).^2));
            % We will only print once by setting 'pp' to 1.
            if pp == 1
                fprintf('Good maximum image contrast = (%.4f / %.4f) / Log sensitivity = (%.4f) \n',marginalContrast, abs(desiredImageContrastCheckCal(index)), log10(1/marginalContrast));
            end
        end
    else
        % Update the order to match the array size. Again, this is not
        % elaborate which will be updated later on.
        pp = pp - nPrimaries;
        
        % PointCloud method
        plot(desiredImageContrastCheckCal,ptCldScreenContrastCheckCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
        plot(desiredImageContrastCheckCal,desiredContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    end
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired image contrast','fontsize',15);
    ylabel('Nominal cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    grid on;
    
    % Add legend.
    switch calculationMethod
        case 'Standard'
            if TESTONEPOINT
                legend('Nominal (Standard)','Desired (Standard)','Cut-off point','location','southeast','fontsize',11);
            else
                legend('Nominal (Standard)','Desired (Standard)', 'location','southeast','fontsize',11);
            end
        case 'PointCloud'
            legend('Nominal (PointCloud)','Desired (PointCloud)','location','southeast','fontsize',11);
        otherwise
    end
end

%% After the meeting (David and Semin, 09/26/23), we added one more plot
if ~exist('desiredImageContrastCheckCal')
    desiredImageContrastCheckCal = theData.spatialGaborTargetContrast * rawMonochromeContrastCheckCal;
end

figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
sgtitle('Desired image contrast vs. Nominal contrast');
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% Again, not fancy loop. Same graph as the above, but different way to look
% at it with updated x-axis as the desired image contrast.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % Save the calculation type per each order.
    calculationMethod = 'Standard';
    
    plot(desiredImageContrastCheckCal,standardScreenContrastMeasuredCheckCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    plot(desiredImageContrastCheckCal,desiredContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    
    % Here we sort the desired contrast (x-axis) in an acsending
    % order so that we can make the line plot correct.
    [desiredImageContrastCheckCalSorted I] = sort(desiredImageContrastCheckCal);
    standardPredictedContrastGaborCalSorted = standardPredictedContrastGaborCal(:,I);
    
    % Set line color.
    lineColorHandle = zeros(1,3);
    lineColorHandle(pp) = 1;
    lineColorAlpha = 0.2;
    lineColorHandle(end+1) = lineColorAlpha;
    plot(desiredImageContrastCheckCalSorted,standardPredictedContrastGaborCalSorted(pp,:),'-','color',lineColorHandle,'linewidth',6);
    
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired image contrast','fontsize',15);
    ylabel('Cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    grid on;
    
    % Add legend.
    legend('Measured (Standard)','Desired (Standard)','Nominal (Standard)', 'location','southeast','fontsize',11);
end
