% SpectralCalCheck
%
% Read in output of SpectralTestCal and make some measurements.
% 
% 10/19/2021  dhb,smo  Started on it
% 10/28/2021  dhb      Add condition name and date to output
% 11/12/2021  dhb      Moving to warmup and measure with steady subprimaries.

%% Initialize.
clear; close all;

%% Parameters
warmupTimeMinutes = 10;
verbose = true;

%% Which condition
%
% This is used to match up with parameters run in SpectralCalCompute
% ['LminusMSmooth' 'ConeIsolating']
conditionName = 'ConeIsolating';

%% Load output of SpectralCalCompute.
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theData = load(testFilename);
end

%% Measure the desired target primaries
%
% This just needs to measure.
%
% Target Spds.
targetPrimarySpd = theData.projectorCalObj.get('P_device');

% Set some variables.
S = theData.S;           % Range of the spectrum.
wls = SToWls(S);         % Wavelength. 
nPrimaries = 3;          % Number of primaries.
nSubprimaries = theData.subprimaryCalObjs{1}.get('nDevices');% Number of subprimaries.
subprimaryNInputLevels = size(theData.subprimaryCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:7 9:15];
nTestPoints = size(theData.thePointCloudContrastCheckCal,2);
T_cones = theData.T_cones;

% Loop and measure all primaries.
%
% IF MEASURE is false, load in the data from a previous run where MEASURE
% was true.
MEASURE = true;
if (MEASURE)
    % Open up projector and radiometer.
    [window,windowRect] = OpenProjectorPlainScreen([1 1 1]');
    OpenSpectroradiometer;
    
    % Set subprimaries to desired value and wait for them to warm up to
    % steady state.
    SetSubprimarySettings(theData.projectorPrimarySettings,'nInputLevels',subprimaryNInputLevels,'projectorMode',true); 
    if (verbose)
        fprintf('Waiting for warmup time of %d minutes ...',warmupTimeMinutes);
    end
    pause(60*warmupTimeMinutes);
    if (verbose)
        fprintf('done.  Measuring.\n');
    end
    
    % Measure
    for pp = 1:nPrimaries
        % Old way. Remove once debugged.
        % isolatingSpdMeasured(:,pp) = MeasureDesiredTargetPrimaries(theData.projectorPrimaryPrimaries(:,pp), ...
        %     theData.subprimaryCalObjs{pp},pp,'projectorMode',true,'measurementOption',true,'verbose',verbose);
        theProjectorOnePrimarySettings = zeros(nPrimaries,1);
        theProjectorOnePrimarySettings(pp) = 1;
        isolatingSpdMeasured(:,pp) = MeasureProjectorPlainScreenSettings(theProjectorOnePrimarySettings,...
            S,window,windowRect,'measurementOption',true,'verbose',verbose);
        clear theProjectorOnePrimarySettings
        
    end
else
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('testImageDataCheck_%s',conditionName),'olderDate',0);
        load(testFilename); 
    else
        error('No file to load');
    end
end

% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in targetPrimarySpd, what we got is in
% isolatingSpdMeasured.
% Plot the spd results.
figure; clf; 
for pp = 1:nPrimaries
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetPrimarySpd(:,pp),'k','LineWidth',3)
    plot(wls,isolatingSpdMeasured(:,pp),'r','LineWidth',2);
    xlabel('Wavelength (nm)');
    ylabel('Spectral power distribution');
    legend('Target','Measured');
    title('Comparison of raw measured and desired spds');
end

%% Get scale factor between target and measured and plot that comparison too
for pp = 1:nPrimaries
    scalePrimaryToTargetFactor(pp) = isolatingSpdMeasured(:,pp)\targetPrimarySpd(:,pp);
    fprintf('Scale factor in measurement for primary %d is %0.3f\n',pp,scalePrimaryToTargetFactor(pp));
end
meanPrimaryScaleFactor = mean(scalePrimaryToTargetFactor);

% Make the plot
figure; clf; 
for pp = 1:nPrimaries
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetPrimarySpd(:,pp),'k','LineWidth',3)
    plot(wls,scalePrimaryToTargetFactor(pp)*isolatingSpdMeasured(:,pp),'r','LineWidth',2);
    xlabel('Wavelength (nm)');
    ylabel('Spectral power distribution');
    legend('Target','Measured');
    title('Comparison of scaled measured and desired spds');
end

%% Set each primary to the settings we loaded in and measure
if (MEASURE)
    projectorCalObj = theData.projectorCalObj;
    theData = rmfield(theData,'projectorCalObj');

    % Set the projector subprimaries here. Now done above.  Remove once
    % debugged.
    % SetSubprimarySettings(theData.projectorPrimarySettings,'nInputLevels',subprimaryNInputLevels,'projectorMode',true);
    
    %% Set the primaries in the calibration to the measured results.
    %
    % It's important to also set the sensor color space, because the
    % transform between sensor/primaries is cached when we set it.
    projectorCalObj.set('P_device',isolatingSpdMeasured);
    SetSensorColorSpace(projectorCalObj,T_cones,S);
    
    %% Optional recompute of target settings
    %
    % We just measured the primaries.  So, we could recompute the contrast
    % check settings based on these measured primaries, rather than the
    % nominal ones we used in SpectralTestCal.  Here is where that happens
    RECOMPUTE = true;
    if (RECOMPUTE)
        % Parameters.
        nQuantizeLevels = theData.nQuantizeLevels;

        % Figure out desired background excitations. The actual background
        % won't be what we computed originally, because the primary spectra
        % aren't the same as their nominal values.  Here we recompute what
        % the background will be when we set the projector to the
        % background settings.  That then lets us compute contrast relative
        % to the background we're going to get. We want to do this from the
        % settings, so that the background isn't mucked up by quantization.
        projectorBgSpd = PrimaryToSpd(projectorCalObj,SettingsToPrimary(projectorCalObj,theData.thePointCloudSettingsCheckCal(:,1)));
        projectorBgExcitations = T_cones * projectorBgSpd;
        
        % Build point cloud
        tic;
        fprintf('Point cloud exhaustive method, setting up cone contrast cloud, this takes a while\n')
        allProjectorIntegersCal = zeros(3,theData.projectorNInputLevels^3);
        idx = 1;
        for ii = 0:(theData.projectorNInputLevels-1)
            for jj = 0:(theData.projectorNInputLevels-1)
                for kk = 0:(theData.projectorNInputLevels-1)
                    allProjectorIntegersCal(:,idx) = [ii jj kk]';
                    idx = idx+1;
                end
            end
        end
        
        % Convert integers to 0-1 reals, quantized
        allProjectorSettingsCal = IntegersToSettings(allProjectorIntegersCal,'nInputLevels',theData.projectorNInputLevels);

        % Get LMS excitations for each triplet of projector settings, and build a
        % point cloud object from these.
        allProjectorExcitations = SettingsToSensor(projectorCalObj,allProjectorSettingsCal);
        allProjectorContrast = ExcitationsToContrast(allProjectorExcitations,projectorBgExcitations);
        allSensorPtCloud = pointCloud(allProjectorContrast');
        
        % Force point cloud setup by finding one nearest neighbor. This is slow,
        % but once it is done subsequent calls are considerably faster.
        findNearestNeighbors(allSensorPtCloud,[0 0 0],1);
        toc
        
        %% Generate some settings values corresponding to known contrasts
        %
        % The reason for this is to measure and check these.  This logic follows
        % how we handled an actual gabor image above. The quantization to
        % nQuantizeLevels isn't strictly needed, but nor is it doing harm.
        rawMonochromeUnquantizedContrastCheckCal = [0 0.25 -0.25 0.5 -0.5 1 -1];
        rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
        theDesiredContrastCheckCal = theData.spatialGaborTargetContrast*theData.targetStimulusContrastDir*rawMonochromeContrastCheckCal;
        theDesiredExcitationsCheckCal = ContrastToExcitation(theDesiredContrastCheckCal,projectorBgExcitations);
        
        % For each check calibration find the settings that
        % come as close as possible to producing the desired excitations.
        %
        % If we measure for a uniform field the spectra corresopnding to each of
        % the settings in the columns of thePointCloudSettingsCheckCal, then
        % compute the cone contrasts with respect to the backgound (0 contrast
        % measurement, first settings), we should approximate the cone contrasts in
        % theDesiredContrastCheckCal.
        fprintf('Point cloud exhaustive method, finding settings\n')
        thePointCloudSettingsCheckCal = zeros(3,size(theDesiredContrastCheckCal,2));
        for ll = 1:size(theDesiredContrastCheckCal,2)
            minIndex = findNearestNeighbors(allSensorPtCloud,theDesiredContrastCheckCal(:,ll)',1);
            thePointCloudSettingsCheckCal(:,ll) = allProjectorSettingsCal(:,minIndex);
        end
        thePointCloudPrimariesCheckCal = SettingsToPrimary(projectorCalObj,thePointCloudSettingsCheckCal);
        thePointCloudSpdCheckCal = PrimaryToSpd(projectorCalObj,thePointCloudPrimariesCheckCal);
        thePointCloudExcitationsCheckCal = SettingsToSensor(projectorCalObj,thePointCloudSettingsCheckCal);
        thePointCloudContrastCheckCal = ExcitationsToContrast(thePointCloudExcitationsCheckCal,projectorBgExcitations);
      
    end
    
    %% Measure contrasts of the settings we computed in SpectralTestCal
    %
    % Measure the contrast points.  We've already got the settings so all we
    % need to do is loop through and set a uniform field to each of the
    % settings in thePointCloudSettingsCheckCal and measure the corresponding
    % spd.
    [thePointCloudSpdMeasured, thePointCloudSettingsIntegers] = MeasureProjectorPlainScreenSettings(thePointCloudSettingsCheckCal,...
        S,window,windowRect,'measurementOption',true,'verbose',verbose);
   
end

%% Make plot of measured versus desired spds.
%
% The desired spds are in thePointCloudSpdCheckCal
%
% Choose to use either raw or scaled measurement spectra.
whichToAnalyze = 'raw';
switch (whichToAnalyze)
    case 'raw'
        thePointCloudSpd = thePointCloudSpdMeasured;
    case 'scaled'
        for tt = 1:nTestPoints
            % Find scale factor to best bring into alignment with target
            scaleSpdToTargetFactor(tt) = thePointCloudSpdMeasured(:,tt)\thePointCloudSpdCheckCal(:,pp);
            thePointCloudSpd(:,tt) = scaleSpdToTargetFactor(tt) * thePointCloudSpdMeasured(:,tt);
        end
    otherwise
        error('Unknown analyze type specified');
end

% Plot it.
figure; clf;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize];
set(gcf,'position',figurePosition);
for tt = 1:nTestPoints
    subplot(round(nTestPoints/2),2,tt); hold on;
    plot(wls,thePointCloudSpdCheckCal(:,tt),'k-','LineWidth',3) % Target spectra
    plot(wls,thePointCloudSpd(:,tt),'r-','LineWidth',2); % Measured spectra
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
    title(sprintf('Test %d %s',tt,whichToAnalyze),'fontsize',16)
end

%% Compute cone contrasts for each spectrum relative to the background
%
% We use the fact that the background settings are in the first column.
%
% 'thePointClousdSpd' will be either raw or scaled spectra based on the
% above 'whichToAnalyze' option.
thePointCloudExcitationsMeasured = T_cones * thePointCloudSpd;
thePointCloudBgExcitationsMeasured = thePointCloudExcitationsMeasured(:,1);
thePointCloudContrastMeasured = ExcitationToContrast(thePointCloudExcitationsMeasured,thePointCloudBgExcitationsMeasured);

% Add the nominal contrasts to be compared in the following graph. These
% should be all lined up on the 45-degree line in the following graph.
addNominalContrast = true;
if (addNominalContrast)
    thePointCloudExcitationsNominal = T_cones * thePointCloudSpdCheckCal;
    thePointCloudBgExcitationsNominal = thePointCloudExcitationsNominal(:,1);
    thePointCloudContrastNominal = ExcitationToContrast(thePointCloudExcitationsNominal,thePointCloudBgExcitationsNominal);
else
end

% Plot measured versus desired contrasts
figure; hold on;
plot(theDesiredContrastCheckCalCal(1,:),thePointCloudContrastMeasured(1,:),'ro','MarkerSize',14,'MarkerFaceColor','r');   % L - measured
plot(theDesiredContrastCheckCal(2,:),thePointCloudContrastMeasured(2,:),'go','MarkerSize',12,'MarkerFaceColor','g');   % M - measured
plot(theDesiredContrastCheckCal(3,:),thePointCloudContrastMeasured(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');   % S - measured
if (addNominalContrast)
    plot(theDesiredContrastCheckCal(1,:),thePointCloudContrastNominal(1,:),'ro','MarkerSize',19);   % L - target
    plot(theDesiredContrastCheckCal(2,:),thePointCloudContrastNominal(2,:),'go','MarkerSize',16);   % M - target
    plot(theDesiredContrastCheckCal(3,:),thePointCloudContrastNominal(3,:),'bo','MarkerSize',14);   % S - target
end
xlabel('Desired contrast');
ylabel('Measured contrast');
axisLim = 0.05;
xlim([-axisLim axisLim]);
ylim([-axisLim axisLim]);
axis('square');
line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
grid on;
if (addNominalContrast)
    legend('L-measured','M-measured','S-measured','L-nominal','M-nominal','S-nominal','location','southeast');
else
    legend('L-measured','M-measured','S-measured','location','southeast');
end
title(sprintf('Desired vs. Measured LMS Contrast, %s',whichToAnalyze));

%% Close projector and save out the measurement data.
if (MEASURE)
    % Close
    CloseProjectorScreen;
    CloseSpectroradiometer;
    
    % Save data with the name containing dayTimestr, so that we can
    % automatically load in the most recent output.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,sprintf('testImageDataCheck_%s_%s',conditionName,dayTimestr));
        save(testFilename,'theData','targetPrimarySpd','isolatingSpdMeasured','projectorCalObj', ...
            'projectorBgSpd','projectorBgExcitations','theDesiredContrastCheckCal','theDesiredExcitationsCheckCal', ...
            'thePointCloudSettingsCheckCal','thePointCloudPrimariesCheckCal','thePointCloudSpdCheckCal','thePointCloudExcitationsCheckCal','thePointCloudContrastCheckCal', ...
            'thePointCloudSettingsIntegers','thePointCloudSpdMeasured','thePointCloudContrastMeasured', ...
            'thePointCloudExcitationsNominal','thePointCloudBgExcitationsNominal','thePointCloudContrastNominal');
    end
end
