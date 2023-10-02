% SpectralCalCheck_ver2
%
% Read in output of SpectralCalCompute and make some measurements.
%
% This is the second version, and it does do calculation of the settings
% with two methods, point cloud and standard methods to compare each other.
%

% History:
%    10/19/2021  dhb,smo  Started on it.
%    10/28/2021  dhb      Add condition name and date to output.
%    11/12/2021  dhb      Moving to warmup and measure with steady subprimaries.
%    10/18/2022  smo      Now it loads the screen primary measurement
%                         results that used for making contrast images so
%                         that we can validate this meausrement results are
%                         fine.
%    09/07/2023  smo      Added an option to use standard method (sensor to
%                         primary).
%    09/25/2023  smo      Cleaned up and commented better.
%    09/26/2023  smo      Added a plot to compare cone contrast vs. image
%                         contrast.

%% Initialize.
clear; close all;

%% Set measurement option and image setting calculation methods.
%
% Set it to 'true' for making new measurements, 'false' for nominal
% calculations. We can choose both for primary measurements
% (MEASUREPRIMARY) and test contrast measurements (MEASURETARGETCONTRAST).
MEASUREPRIMARY = false;
MEASURETARGETCONTRAST = false;

% We can calculate the image settings with both PointCloud and Standard
% methods for comparison.
POINTCLOUD = true;
STANDARD = true;

% Control the messages and plots.
verbose = true;

%% Load the image data.
%
% Set the image type to load.
%
% LOADIMAGETYPE to 'new' will load the image data from the results made
% from the code SpectralCalCompute.m.
%
% LOADIMAGETYPE to 'experiment' will load the image data from the SACCSFA
% validation file. This option has been newly added for this new analysis
% comparing between PointCloud and Standard methods.
LOADIMAGETYPE = 'experiment';

% Load the image data here.
switch LOADIMAGETYPE
    % Load new image data.
    case 'new'
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = getpref('SpatioSpectralStimulator','SACCData');
            testFilename = GetMostRecentFileName(fullfile(testFiledir,'CheckCalibration'),...
                'testImageData');
            theData = load(testFilename);
        end
        
        % Load the SACCSFA validation image data.
    case 'experiment'
        % Set which file to load.
        %
        % Default option is to load the most recent validation file
        % ('olderDate' = 0). Otherwise, set it as you want. For example, if
        % olderDate is set to 5, it will search the fifth file to the most
        % recent one.
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
end

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

if isfield(theData,{'spatialGaborTargetContrast','targetScreenPrimaryContrast','targetLambda'})
    spatialGaborTargetContrast = theData.spatialGaborTargetContrast;
    targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
    targetLambda = theData.targetLambda;
end

%% Open up projector and spectroradiometer.
%
% This happens only either when measuring primaries or test contrasts.
if or(MEASURETARGETCONTRAST,MEASUREPRIMARY)
    % Open the projector.
    initialScreenSettings = [1 1 1]';
    [window,windowRect] = OpenPlainScreen(initialScreenSettings);
    % Connect to spectroradiometer.
    OpenSpectroradiometer;
    % Set subprimary settings
    SetChannelSettings(theData.screenPrimarySettings,'nInputLevels',channelNInputLevels);
end

%% Measure primaries here. We can load it too if there is a saved file.
if (MEASUREPRIMARY)
    % Measure.
    for pp = 1:nPrimaries
        theScreenOnePrimarySettings = zeros(nPrimaries,1);
        theScreenOnePrimarySettings(pp) = 1;
        targetScreenSpdMeasured(:,pp) = MeasurePlainScreenSettings(theScreenOnePrimarySettings,...
            S,window,windowRect,'measurementOption',true,'verbose',verbose);
        clear theScreenOnePrimarySettings;
    end
    primaryFilename = 'none';
    
    % Load the measurement results.
elseif (~MEASUREPRIMARY)
    % If we load the data from the experiment, we already saved the
    % measured primary in the file, so we skip the part loading the
    % measured primary file.
    if strcmp(LOADIMAGETYPE,'experiment')
        targetScreenSpdMeasured = theData.screenCalObj.cal.P_device;
    else
        % Otherwise, load the measured primary file.
        if (ispref('SpatioSpectralStimulator','SACCData'))
            % Load different file name according to 'normal' set or 'high' test
            % image contrast sets. We set the primary contrast as 0.07 for
            % normal image, and 0.10 for high image.
            targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
            if (targetScreenPrimaryContrast == 0.10)
                primaryContrast = 'high';
            else
                primaryContrast = 'normal';
            end
            
            % File name convention to load the file.
            switch primaryContrast
                case 'normal'
                    filenamePart = 'targetScreenSpdMeasured_2';
                case 'high'
                    filenamePart = 'targetScreenSpdMeasured_high';
            end
            
            % Load file.
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages','MeasurementData');
            primaryFilename = GetMostRecentFileName(testFiledir,append(filenamePart),'olderdate',olderDate);
            load(primaryFilename);
            fprintf('Measurement file found, so skipping primary measurement! \n');
        else
            error('No file to load');
        end
    end
end

%% Compare what we wanted for primaries versus what we got.
%
% What we want is in 'targetScreenSpd', what we got is in
% 'targetScreenSpdMeasured'.
%
% This part will be run only when we measure primary.
if (MEASUREPRIMARY)
    figure; clf;
    for pp = 1:nPrimaries
        subplot(nPrimaries,1,pp); hold on;
        plot(wls,targetScreenSpd(:,pp),'k','LineWidth',3)
        plot(wls,targetScreenSpdMeasured(:,pp),'r','LineWidth',2);
        xlabel('Wavelength (nm)');
        ylabel('Spectral power distribution');
        legend('Target','Measured');
        title('Comparison of raw measured and desired spds');
    end
    
    % Get scale factor between target and measured and plot that comparison too
    for pp = 1:nPrimaries
        scalePrimaryToTargetFactor(pp) = targetScreenSpdMeasured(:,pp)\targetScreenSpd(:,pp);
        fprintf('Scale factor in measurement for primary %d is %0.3f\n',pp,scalePrimaryToTargetFactor(pp));
    end
    meanPrimaryScaleFactor = mean(scalePrimaryToTargetFactor);
    
    % Make the plot
    figure; clf;
    for pp = 1:nPrimaries
        subplot(nPrimaries,1,pp); hold on;
        plot(wls,targetScreenSpd(:,pp),'k','LineWidth',3)
        plot(wls,scalePrimaryToTargetFactor(pp)*targetScreenSpdMeasured(:,pp),'r','LineWidth',2);
        xlabel('Wavelength (nm)');
        ylabel('Spectral power distribution');
        legend('Target','Measured');
        title('Comparison of scaled measured and desired spds');
    end
end

%% Set each primary to the settings we loaded in and measure
screenCalObj = theData.screenCalObj;
theData = rmfield(theData,'screenCalObj');

%% Set the primaries in the calibration to the measured results.
%
% It's important to also set the sensor color space, because the
% transform between sensor/primaries is cached when we set it.
screenCalObj.set('P_device',targetScreenSpdMeasured);
SetSensorColorSpace(screenCalObj,T_cones,S);

%% Get target settings.
%
% This part is the same for both point cloud and standard methods as it is
% getting 'target' settings.
%
% We just measured the primaries. So, we could recompute the contrast
% check settings based on these measured primaries, rather than the
% nominal ones we used in SpectralTestCal. Here is where that happens.
%
% Generate some settings values corresponding to known contrasts.
%
% The reason for this is to measure and check these. This logic follows
% how we handled an actual gabor image above. The quantization to
% nQuantizeLevels isn't strictly needed, but nor is it doing harm.
nQuantizeLevels = theData.nQuantizeLevels;

% Figure out desired background excitations. The actual background
% won't be what we computed originally, because the primary spectra
% aren't the same as their nominal values. Here we recompute what
% the background will be when we set the screen to the
% background settings. That then lets us compute contrast relative
% to the background we're going to get. We want to do this from the
% settings, so that the background isn't mucked up by quantization.
screenBgSpd = PrimaryToSpd(screenCalObj,SettingsToPrimary(screenCalObj,theData.ptCldScreenSettingsCheckCal(:,1)));
screenBgExcitations = T_cones * screenBgSpd;

% Set test contrasts. We used to use a total of 15 test contrats, but when
% we do nominal settings comparison, we set far more test contrasts.
if (MEASURETARGETCONTRAST)
    % These are 15 contrasts that we used for test image vaildation.
    if isfield(theData,'rawMonochromeUnquantizedContrastCheckCal')
        rawMonochromeUnquantizedContrastCheckCal = theData.rawMonochromeUnquantizedContrastCheckCal;
    else
        rawMonochromeUnquantizedContrastCheckCal = [0 0.05 -0.05 0.10 -0.10 0.15 -0.15 0.20 -0.20 0.25 -0.25 0.5 -0.5 1 -1];
    end
    
elseif (~MEASURETARGETCONTRAST)
    % When we use only nominal comparison skipping the measurements, we set
    % the contrasts a lot more. This setting generates 202 test contrasts.
    % It is fine enough to decide the maximum contrast that generates
    % 'good' image.
    contrastPos = [0:0.01:1];
    contrastNeg = [0:-0.01:-1];
    rawMonochromeUnquantizedContrastCheckCal = [contrastPos contrastNeg];
end

% Calculate the desired cone cone contrasts.
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
desiredImageContrastCheckCal = theData.spatialGaborTargetContrast * rawMonochromeContrastCheckCal;
desiredContrastCheckCal = theData.spatialGaborTargetContrast * theData.targetStimulusContrastDir * rawMonochromeContrastCheckCal;
desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,screenBgExcitations);

%% Here we calculate the settings to achieve the desired contrast.
%
% 1) Point cloud method.
if (POINTCLOUD)
    % Optional recompute of target settings.
    RECOMPUTE = true;
    
    if (RECOMPUTE)
        % Set up point cloud for finding best settings.
        [contrastPtCld,ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',verbose);
        
        % For each check calibration find the settings that come as close
        % as possible to producing the desired excitations.
        %
        % If we measure for a uniform field the spectra corresopnding to
        % each of the settings in the columns of
        % ptCldScreenSettingsCheckCal, then compute the cone contrasts with
        % respect to the backgound (0 contrast measurement, first
        % settings), we should approximate the cone contrasts in
        % 'desiredContrastCheckCal'.
        fprintf('Point cloud exhaustive method, finding settings\n')
        ptCldScreenSettingsCheckCal = SettingsFromPointCloud(contrastPtCld,desiredContrastCheckCal,ptCldSettingsCal);
        ptCldScreenPrimariesCheckCal = SettingsToPrimary(screenCalObj,ptCldScreenSettingsCheckCal);
        ptCldScreenSpdCheckCal = PrimaryToSpd(screenCalObj,ptCldScreenPrimariesCheckCal);
        ptCldScreenExcitationsCheckCal = SettingsToSensor(screenCalObj,ptCldScreenSettingsCheckCal);
        ptCldScreenContrastCheckCal = ExcitationsToContrast(ptCldScreenExcitationsCheckCal,screenBgExcitations);
    end
end

% 2) Standard method.
if (STANDARD)
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
end

%% Compare the nominal contrasts to compate Standard method to PointCloud.
%
% Get Point cloud settings when it is not recomputed. As we increased the
% number of test contrasts when we compare the nominal contrasts, so
% RECOMPUTE should be set as 'true' when we are testing more than 15
% contrats, which were used in validation for SACCSFA experiment.
if (~RECOMPUTE)
    ptCldScreenContrastCheckCal = theData.ptCldScreenContrastCheckCal;
end

%% Plot 1) Desired vs. Nominal contrasts.
figure; hold on;
figurePosition = [0 0 1300 700];
set(gcf,'position',figurePosition);
sgtitle('Nominal contrast: Standard vs. Point cloud methods');
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% Make a loop for plotting each cone case. We decided to plot the results
% separately over the calculation methods, Standard and PointCloud. Making
% a loop in this way is not the most elaborate, so we may want to update
% this part later on.
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
    switch calculationMethod
        case 'Standard'
            % Standard method.
            plot(desiredContrastCheckCal(pp,:),standardPredictedContrastGaborCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
        case 'PointCloud'
            % Update the order to match the array size. Again, this is not
            % elaborate which will be updated later on.
            pp = pp - nPrimaries;
            % PointCloud method
            plot(desiredContrastCheckCal(pp,:),ptCldScreenContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    end
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired cone contrast','fontsize',15);
    ylabel('Nominal cone contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    
    % Add legend.
    switch calculationMethod
        case 'Standard'
            legend('Standard','location','southeast','fontsize',13);
        case 'PointCloud'
            legend('PointCloud','location','southeast','fontsize',13);
        otherwise
    end
end

% From the above figure, we visually searched that the number of index
% (183) of the test contrasts was about the marginal contrast that makes a
% good test images. Its cone contrast was [L M S] = [-0.0427 0.0369
% -0.0028] which generates the image contrast of 0.0565. The other side
% (index=101) was good [L M S] = [0.0493 -0.0496 -0.0003], which generates
% the image contrast of 0.0699. This is for normal image set.
%
% Likewise, for high image set, the good test image index was (185) and its
% cone contrast was [L M S] = [-0.0579 0.0590 -0.0003] which generates the
% image contrast of 0.0827. Also, the other side was good (index = 101),
% where its cone contrast was [L M S] = [0.0694 -0.0717 -0.0012], its test
% image contrast is 0.0998.

%% Plot 2) Nominal cone contrast over desired test image contrast.
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
        
        % Testing one point. Change the index from 1 to 202 (number of
        % nominal contrasts to test) so that you can visually check from
        % when the contrast reproduction goes not great.
        %
        % For now, we set the cone contrast as the same index within the
        % test contrats.
        TESTONEPOINT = true;
        if(TESTONEPOINT)
            index = 185;
            plot(desiredImageContrastCheckCal(index),standardPredictedContrastGaborCal(pp,index),'o','MarkerSize',14,'MarkerFaceColor','c');
            fprintf('Good maximum image contrast = (%.4f) / Log sensitivity = (%.4f) \n',desiredImageContrastCheckCal(index), log10(1/desiredImageContrastCheckCal(index)));
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
            legend('Nominal (Standard)','Desired (Standard)', 'location','southeast','fontsize',11);
        case 'PointCloud'
            legend('Nominal (PointCloud)','Desired (PointCloud)','location','southeast','fontsize',11);
        otherwise
    end
end

%% Plot 3) Direct comparison of nominal contrasts between Standard vs.
% PointCloud.
figure; clf;
set(gcf,'position',figurePosition);
sgtitle('Direct comparison of nominal contrasts: Standard vs. Point cloud methods');

% Here we make a loop again for each cone case.
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    plot(standardPredictedContrastGaborCal(pp,:),ptCldScreenContrastCheckCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    title(titleHandles{pp},'fontsize',15);
    xlabel('Nominal contrast (Standard)','fontsize',15);
    ylabel('Nominal contrast (PointCloud)','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
end

%% From here, we measure the test contrasts using spectroradiometer.
%
% Measure the contrast points. We've already got the settings so all we
% need to do is loop through and set a uniform field to each of the
% settings in ptCldScreenSettingsCheckCal and measure the corresponding
% spd.
%
% We updated this code to measure the test contrasts from two different
% methods, PointCloud and Standard at the same time. As we set the test
% contrasts as 15, so here we will make total 30 measurements (15 contrasts
% x 2 different methods).
if (MEASURETARGETCONTRAST)
    % Settings using Point cloud method.
    if (POINTCLOUD)
        disp('Measurements wil begin! - (Settings: Point cloud)');
        [ptCldScreenSpdMeasuredCheckCal, ptCldScreenSettingsIntegersCheckCal] = MeasurePlainScreenSettings(ptCldScreenSettingsCheckCal,...
            S,window,windowRect,'measurementOption',true,'verbose',verbose);
        disp('Measurements finished - (Settings: Point cloud)');
    end
    
    % Settings using the standard method.
    if (STANDARD)
        disp('Measurements wil begin! - (Settings: Standard)');
        [standardScreenSpdMeasuredCheckCal, standardScreenSettingsIntegersCheckCal] = MeasurePlainScreenSettings(standardSettingsGaborCal,...
            S,window,windowRect,'measurementOption',true,'verbose',verbose);
        disp('Measurements finished - (Settings: Standard)');
    end
end

%% Make plot of measured versus desired spds.
%
% The desired spds from  are in ptCldScreenSpdCheckCal. The measured spds
% from PointCloud is in ptCldSpd and the measured spds from Standard method
% is in standardSpd. We just shortened the variable names here not to be
% confused.
if (MEASURETARGETCONTRAST)
    ptCldSpd = ptCldScreenSpdMeasuredCheckCal;
    standardSpd = standardScreenSpdMeasuredCheckCal;
    
    % Plot it.
    figure; clf;
    figureSize = 1000;
    figurePosition = [1200 300 figureSize figureSize];
    set(gcf,'position',figurePosition);
    
    nTestPoints = size(rawMonochromeUnquantizedContrastCheckCal,2);
    for tt = 1:nTestPoints
        subplot(round(nTestPoints/2),2,tt); hold on;
        % Target spectra.
        plot(wls,ptCldScreenSpdCheckCal(:,tt),'k-','LineWidth',3);
        % Measured spectra (PointCloud).
        plot(wls,ptCldSpd(:,tt),'r-','LineWidth',2);
        % Measure spectra (Standard).
        plot(wls,standardSpd(:,tt),'b--','LineWidth',2);
        xlabel('Wavelength (nm)')
        ylabel('Spectral power distribution')
        legend('Target','Measured(PT)','Measured(ST)')
        title(sprintf('Test %d',tt),'fontsize',16)
    end
end

%% Compute cone contrasts for each spectrum relative to the background.
%
% This part only runs when we measured the test contrasts.
if (MEASURETARGETCONTRAST)
    % 1) Point cloud
    % We use the fact that the background settings are in the first column.
    ptCldExcitationsMeasured = T_cones * ptCldSpd;
    ptCldBgExcitationsMeasured = ptCldExcitationsMeasured(:,1);
    ptCldScreenContrastMeasuredCheckCal = ExcitationToContrast(ptCldExcitationsMeasured,ptCldBgExcitationsMeasured);
    
    % Add the nominal contrasts to be compared in the following graph. These
    % should be all lined up on the 45-degree line in the following graph.
    ptCldExcitationsNominal = T_cones * ptCldScreenSpdCheckCal;
    ptCldBgExcitationsNominal = ptCldExcitationsNominal(:,1);
    ptCldContrastNominal = ExcitationToContrast(ptCldExcitationsNominal,ptCldBgExcitationsNominal);
    
    % Plot measured versus desired contrasts
    figure; hold on;
    % Measured contrast.
    plot(desiredContrastCheckCal(1,:),ptCldScreenContrastMeasuredCheckCal(1,:),'ro','MarkerSize',14,'MarkerFaceColor','r');   % L - measured
    plot(desiredContrastCheckCal(2,:),ptCldScreenContrastMeasuredCheckCal(2,:),'go','MarkerSize',12,'MarkerFaceColor','g');   % M - measured
    plot(desiredContrastCheckCal(3,:),ptCldScreenContrastMeasuredCheckCal(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');   % S - measured
    
    % Nominal contrast.
    plot(desiredContrastCheckCal(1,:),ptCldContrastNominal(1,:),'ro','MarkerSize',19);   % L - target
    plot(desiredContrastCheckCal(2,:),ptCldContrastNominal(2,:),'go','MarkerSize',16);   % M - target
    plot(desiredContrastCheckCal(3,:),ptCldContrastNominal(3,:),'bo','MarkerSize',14);   % S - target
    
    xlabel('Desired contrast','fontsize',15);
    ylabel('Measured contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('L-measured','M-measured','S-measured','L-nominal','M-nominal','S-nominal','location','southeast','fontsize',12);
    title(sprintf('Desired vs. Measured LMS Contrast, %s','Point cloud'));
    
    % 2) Standard method
    % We use the fact that the background settings are in the first column.
    standardExcitationsMeasured = T_cones * standardSpd;
    standardBgExcitationsMeasured = standardExcitationsMeasured(:,1);
    standardScreenContrastMeasuredCheckCal = ExcitationToContrast(standardExcitationsMeasured,standardBgExcitationsMeasured);
    
    % Add the nominal contrasts to be compared in the following graph.
    standardContrastNominal = standardPredictedContrastGaborCal;
    
    % Plot measured versus desired contrasts
    figure; hold on;
    % Measured contrast.
    plot(desiredContrastCheckCal(1,:),standardScreenContrastMeasuredCheckCal(1,:),'ro','MarkerSize',14,'MarkerFaceColor','r');   % L - measured
    plot(desiredContrastCheckCal(2,:),standardScreenContrastMeasuredCheckCal(2,:),'go','MarkerSize',12,'MarkerFaceColor','g');   % M - measured
    plot(desiredContrastCheckCal(3,:),standardScreenContrastMeasuredCheckCal(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');   % S - measured
    
    % Nominal contrast.
    plot(desiredContrastCheckCal(1,:),standardContrastNominal(1,:),'ro','MarkerSize',19);   % L - target
    plot(desiredContrastCheckCal(2,:),standardContrastNominal(2,:),'go','MarkerSize',16);   % M - target
    plot(desiredContrastCheckCal(3,:),standardContrastNominal(3,:),'bo','MarkerSize',14);   % S - target
    
    xlabel('Desired contrast','fontsize',15);
    ylabel('Measured contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('L-measured','M-measured','S-measured','L-nominal','M-nominal','S-nominal','location','southeast','fontsize',12);
    title(sprintf('Desired vs. Measured LMS Contrast, %s','Standard method'));
end

%% Close screen and save out the measurement data.
%
% If we measured new target contrasts, we close the projector and
% spectroradiometer, and save the results.
if (MEASURETARGETCONTRAST)
    % Close the projector and spectroradiometer.
    CloseScreen;
    CloseSpectroradiometer;
    
    % Save data with the name containing dayTimestr, so that we can
    % automatically load in the most recent output.
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,'CheckCalibration',sprintf('testImageDataCheck_postExp_%s',dayTimestr));
        save(testFilename,'theData','targetScreenSpd','targetScreenSpdMeasured','screenCalObj', ...
            'screenBgSpd','screenBgExcitations','desiredContrastCheckCal','desiredExcitationsCheckCal', ...
            'ptCldScreenSettingsCheckCal','ptCldScreenPrimariesCheckCal','ptCldScreenSpdCheckCal','ptCldScreenExcitationsCheckCal','ptCldScreenContrastCheckCal', ...
            'ptCldScreenSettingsIntegersCheckCal','ptCldScreenSpdMeasuredCheckCal','ptCldScreenContrastMeasuredCheckCal', ...
            'ptCldExcitationsNominal','ptCldBgExcitationsNominal','ptCldContrastNominal',...
            'primaryFilename','targetScreenPrimaryContrast','targetLambda','spatialGaborTargetContrast',...
            'standardScreenSpdMeasuredCheckCal','standardScreenContrastMeasuredCheckCal','standardDesiredSpdGaborCal','standardContrastNominal',...
            'standardDesiredSpdGaborCal','standardExcitationsMeasured');
    end
    disp('Data has been saved successfully!');
end

%% This part is from SpectralCalAnalyze.
%
% This shows basically the same results, but we put it here, which are
% origianlly from the analyzing code, SpectralCalAnalyze.m.
if (MEASURETARGETCONTRAST)
    % 1) Point cloud.
    %
    % Set figure size and position.
    contrastFig = figure; hold on;
    figureSize = 1000;
    figurePosition = [1200 300 figureSize figureSize/3];
    set(gcf,'position',figurePosition);
    
    % Plot measured versus desired contrasts.
    axisLim = 0.10;
    theColors = ['r' 'g' 'b'];
    sgtitle('Point cloud method');
    for pp = 1:nPrimaries
        subplot(1,nPrimaries,pp); hold on;
        plot(desiredContrastCheckCal(pp,:),ptCldScreenContrastMeasuredCheckCal(pp,:),[theColors(pp) 'o'],'MarkerSize',14,'MarkerFaceColor',theColors(pp));
        plot(desiredContrastCheckCal(pp,:),ptCldContrastNominal(pp,:), [theColors(pp) 'o'],'MarkerSize',17);
        plot(desiredContrastCheckCal(pp,1),ptCldScreenContrastMeasuredCheckCal(pp,1),'ko','MarkerSize',14,'MarkerFaceColor','k');
        plot(desiredContrastCheckCal(pp,1),ptCldContrastNominal(pp,1), 'ko','MarkerSize',17);
        
        plot([-1 1],[-1 1],'k');
        xlim([-axisLim axisLim]);
        ylim([-axisLim axisLim]);
        axis('square');
        xlabel('Desired contrast','fontsize',15);
        ylabel('Measured contrast','fontsize',15);
        legend({'Measured','Nominal'},'location','southeast');
        title(sprintf('Cone class %d',pp),'fontsize',15);
        grid on;
    end
    
    % 2) Standard method.
    %
    % Set figure size and position.
    contrastFig = figure; hold on;
    figureSize = 1000;
    figurePosition = [1200 300 figureSize figureSize/3];
    set(gcf,'position',figurePosition);
    
    % Plot measured versus desired contrasts.
    axisLim = 0.10;
    theColors = ['r' 'g' 'b'];
    sgtitle('Standard method');
    for pp = 1:nPrimaries
        subplot(1,nPrimaries,pp); hold on;
        plot(desiredContrastCheckCal(pp,:),standardScreenContrastMeasuredCheckCal(pp,:),[theColors(pp) 'o'],'MarkerSize',14,'MarkerFaceColor',theColors(pp));
        plot(desiredContrastCheckCal(pp,:),standardContrastNominal(pp,:), [theColors(pp) 'o'],'MarkerSize',17);
        plot(desiredContrastCheckCal(pp,1),standardScreenContrastMeasuredCheckCal(pp,1),'ko','MarkerSize',14,'MarkerFaceColor','k');
        plot(desiredContrastCheckCal(pp,1),standardContrastNominal(pp,1), 'ko','MarkerSize',17);
        
        plot([-1 1],[-1 1],'k');
        xlim([-axisLim axisLim]);
        ylim([-axisLim axisLim]);
        axis('square');
        xlabel('Desired contrast','fontsize',15);
        ylabel('Measured contrast','fontsize',15);
        legend({'Measured','Nominal'},'location','southeast');
        title(sprintf('Cone class %d',pp),'fontsize',15);
        grid on;
    end
    
    %% After the meeting (David and Semin, 09/26/23), we added one more plot
    % comparing the desired image contrast vs. cone contrasts.
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
            plot(desiredImageContrastCheckCal,standardScreenContrastMeasuredCheckCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
            plot(desiredImageContrastCheckCal,standardPredictedContrastGaborCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
        else
            % Update the order to match the array size. Again, this is not
            % elaborate which will be updated later on.
            pp = pp - nPrimaries;
            
            % PointCloud method
            plot(desiredImageContrastCheckCal,ptCldScreenContrastMeasuredCheckCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
            plot(desiredImageContrastCheckCal,ptCldScreenContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
        end
        title(titleHandles{pp},'fontsize',15);
        xlabel('Desired image contrast','fontsize',15);
        ylabel('Cone contrast','fontsize',15);
        axisLim = 0.10;
        xlim([-axisLim axisLim]);
        ylim([-axisLim axisLim]);
        axis('square');
        grid on;
        
        % Add legend.
        switch calculationMethod
            case 'Standard'
                legend('Measured (Standard)','Nominal (Standard)', 'location','southeast','fontsize',11);
            case 'PointCloud'
                legend('Measured (PointCloud)','Nominal (PointCloud)','location','southeast','fontsize',11);
            otherwise
        end
    end
end
