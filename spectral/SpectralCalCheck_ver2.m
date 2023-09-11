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

%% Initialize.
clear; close all;

%% Set measurement and setting calculation options.
verbose = true;
MEASUREPRIMARY = false;
MEASURETARGETCONTRAST = false;
POINTCLOUD = true;
STANDARD = true;

%% Load the image data.
%
% You can load image either from fresh test image or from the available image
% from SACCSFA experiment.
LOADIMAGETYPE = 'experiment';
switch LOADIMAGETYPE
    case 'new'
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = getpref('SpatioSpectralStimulator','SACCData');
            testFilename = GetMostRecentFileName(fullfile(testFiledir,'CheckCalibration'),...
                'testImageData');
            theData = load(testFilename);
        end
        
    case 'experiment'
        % Read out the calibration data.
        olderDate = 0;
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'CheckCalibration');
            testFilename = GetMostRecentFileName(testFiledir,'testImageDataCheck_','olderDate',olderDate+1);
            theData = load(testFilename);
        end
        % Match the data format to run it smoothly.
        screenCalObj_temp = theData.screenCalObj;
        theData = theData.theData;
        theData.screenCalObj = screenCalObj_temp;
        
        % Get the date of the experiment.
        numExtract = regexp(testFilename,'\d+','match');
        strYear = numExtract{1};
        strMonth = numExtract{2};
        strDay = numExtract{3};
        measureDate = sprintf('_%s-%s-%s',strYear,strMonth,strDay);
end

%% Set variables here.
targetScreenSpd = theData.screenCalObj.get('P_device');
S = theData.S;
wls = SToWls(S);
nPrimaries = 3;
nChannels = theData.channelCalObjs{1}.get('nDevices');
channelNInputLevels = size(theData.channelCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:15];
nTestPoints = size(theData.ptCldScreenContrastCheckCal,2);
T_cones = theData.T_cones;

% Newly added variables.
if isfield(theData,{'spatialGaborTargetContrast','targetScreenPrimaryContrast','targetLambda'})
    spatialGaborTargetContrast = theData.spatialGaborTargetContrast;
    targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
    targetLambda = theData.targetLambda;
end

%% Open up projector and spectroradiometer.
if or(MEASURETARGETCONTRAST,MEASUREPRIMARY)
    % Open the projector.
    initialScreenSettings = [1 1 1]';
    [window,windowRect] = OpenPlainScreen(initialScreenSettings);
    % Connect to spectroradiometer.
    OpenSpectroradiometer;
    % Set subprimary settings
    SetChannelSettings(theData.screenPrimarySettings,'nInputLevels',channelNInputLevels);
end

%% Measure primaries here. We can load it too if there is a file saved.
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
            % image contrast sets.
            targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
            if (targetScreenPrimaryContrast > 0.07)
                primaryContrast = 'high';
            else
                primaryContrast = 'normal';
            end
            
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

% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in 'targetScreenSpd', what we got is in
% 'targetScreenSpdMeasured'.
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

%% Get scale factor between target and measured and plot that comparison too
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

if isfield(theData,'rawMonochromeUnquantizedContrastCheckCal')
    rawMonochromeUnquantizedContrastCheckCal = theData.rawMonochromeUnquantizedContrastCheckCal;
else
    rawMonochromeUnquantizedContrastCheckCal = [0 0.05 -0.05 0.10 -0.10 0.15 -0.15 0.20 -0.20 0.25 -0.25 0.5 -0.5 1 -1];
end

% We set more contrasts to test if we skip the measurements.
if ~and(MEASUREPRIMARY,MEASURETARGETCONTRAST)
    contrastPos = [0:0.01:1];
    contrastNeg = [0:-0.01:-1];
    rawMonochromeUnquantizedContrastCheckCal = [contrastPos contrastNeg];
end

% Calculate the desired cone cone contrasts.
rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
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

% 2) Standard method
if (STANDARD)
    % Get primaries using standard calibration code, and desired spd without
    % quantizing.
    standardPrimariesGaborCal = SensorToPrimary(screenCalObj,desiredExcitationsCheckCal);
    standardDesiredSpdGaborCal = PrimaryToSpd(screenCalObj,standardPrimariesGaborCal);
    
    % Gamma correct and quantize (if gamma method set to 2 above; with gamma
    % method set to zero there is no quantization).  Then convert back from
    % the gamma corrected settings.
    standardSettingsGaborCal = PrimaryToSettings(screenCalObj,standardPrimariesGaborCal);
    standardPredictedPrimariesGaborCal = SettingsToPrimary(screenCalObj,standardSettingsGaborCal);
    standardPredictedExcitationsGaborCal = PrimaryToSensor(screenCalObj,standardPredictedPrimariesGaborCal);
    standardPredictedContrastGaborCal = ExcitationsToContrast(standardPredictedExcitationsGaborCal,screenBgExcitations);
end

%% Compare the nominal contrasts between standard and point cloud methods.
%
% Get Point cloud settings unless it is re-computed.
if (~RECOMPUTE)
    ptCldScreenContrastCheckCal = theData.ptCldScreenContrastCheckCal;
end

% Plot it here.
figure; hold on;
figurePosition = [0 0 1500 500];
set(gcf,'position',figurePosition);
sgtitle('Nominal contrast: Standard vs. Point cloud methods');
titleHandles = {'L-cone', 'M-cone', 'S-cone'};
markerColorHandles = {'r','g','b'};

% testIndex = 183 for normal image set when we set the interval of 0.01.
% The highest cone cone trast [L M S] = [-0.0427 0.0369 -0.0028] -> image
% contrast = 0.0565.
% testIndex = 101 for normal image set. The other side of good. [0.0493
% -0.0496 -0.0003] -> 0.0699
%
% textIndex = 185 for high image set.
% The highest cone cone trast [L M S] = [-0.0579 0.0590 -0.0003] -> image
% contrast = 0.0827
% testIndex = 101 for high image set. The other side of good. [0.0694
% -0.0717 -0.0012] -> 0.0998
testIndex = 101;
for pp = 1:nPrimaries
    subplot(1,3,pp); hold on;
    
    % Standard method.
    plot(desiredContrastCheckCal(pp,:),standardPredictedContrastGaborCal(pp,:),'o','MarkerSize',14,'MarkerFaceColor',markerColorHandles{pp});
    
    % PointCloud method
    plot(desiredContrastCheckCal(pp,:),ptCldScreenContrastCheckCal(pp,:),'o','MarkerSize',17,'MarkerEdgeColor',markerColorHandles{pp});
    
    % Test - Standard method.
    plot(desiredContrastCheckCal(pp,1:testIndex),standardPredictedContrastGaborCal(pp,1:testIndex),'o','MarkerSize',14,'MarkerFaceColor','c'); 
    
    title(titleHandles{pp},'fontsize',15);
    xlabel('Desired contrast','fontsize',15);
    ylabel('Nominal contrast','fontsize',15);
    axisLim = 0.10;
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
    grid on;
    legend('Standard','PointCloud','location','southeast','fontsize',13);
end

% Another way to see it.
figure; clf;
set(gcf,'position',figurePosition);
sgtitle('Direct comparison of nominal contrasts: Standard vs. Point cloud methods');
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

%% Measure contrasts of the settings we computed in SpectralTestCal.
%
% Measure the contrast points. We've already got the settings so all we
% need to do is loop through and set a uniform field to each of the
% settings in ptCldScreenSettingsCheckCal and measure the corresponding
% spd.
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
% The desired spds are in ptCldScreenSpdCheckCal
if (MEASURETARGETCONTRAST)
    ptCldSpd = ptCldScreenSpdMeasuredCheckCal;
    standardSpd = standardScreenSpdMeasuredCheckCal;
    
    % Plot it.
    figure; clf;
    figureSize = 1000;
    figurePosition = [1200 300 figureSize figureSize];
    set(gcf,'position',figurePosition);
    for tt = 1:nTestPoints
        subplot(round(nTestPoints/2),2,tt); hold on;
        plot(wls,ptCldScreenSpdCheckCal(:,tt),'k-','LineWidth',3) % Target spectra
        plot(wls,ptCldSpd(:,tt),'r-','LineWidth',2); % Measured spectra - point cloud
        plot(wls,standardSpd(:,tt),'b--','LineWidth',2); % Measured spectra - standard
        xlabel('Wavelength (nm)')
        ylabel('Spectral power distribution')
        legend('Target','Measured(PT)','Measured(ST)')
        title(sprintf('Test %d',tt),'fontsize',16)
    end
end

%% Compute cone contrasts for each spectrum relative to the background
%
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
if (MEASURETARGETCONTRAST)
    % Close.
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
if (MEASURETARGETCONTRAST)
    % 1) Point cloud
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
    
    % 2) Standard method
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
end