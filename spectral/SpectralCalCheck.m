% SpectralCalCheck
%
% Read in output of SpectralTestCal and make some measurements.
% 

% History:
%    10/19/2021  dhb,smo  Started on it
%    10/28/2021  dhb      Add condition name and date to output
%    11/12/2021  dhb      Moving to warmup and measure with steady subprimaries.

%% Initialize.
clear; close all;

%% Parameters
warmupTimeMinutes = 0;
verbose = true;
MEASURE = true;

%% Which condition
%
% This is used to match up with parameters run in SpectralCalCompute
% ['LminusMSmooth' 'ConeIsolating']
conditionName = 'LminusMSmooth';

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
targetScreenSpd = theData.screenCalObj.get('P_device');

% Set some variables.
S = theData.S;           % Range of the spectrum.
wls = SToWls(S);         % Wavelength. 
nPrimaries = 3;          % Number of primaries.
nChannels = theData.channelCalObjs{1}.get('nDevices');% Number of subprimaries.
channelNInputLevels = size(theData.channelCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:7 9:15];
nTestPoints = size(theData.ptCldScreenContrastCheckCal,2);
T_cones = theData.T_cones;

% Loop and measure all primaries.
%
% IF MEASURE is false, load in the data from a previous run where MEASURE
% was true.
if (MEASURE)
    % Open up screen and radiometer.
    [window,windowRect] = OpenPlainScreen([1 1 1]');
    OpenSpectroradiometer;
    
    % Set subprimaries to desired value and wait for them to warm up to
    % steady state.
    SetChannelSettings(theData.screenPrimarySettings,'nInputLevels',channelNInputLevels); 
    if (verbose)
        fprintf('Waiting for warmup time of %d minutes ...',warmupTimeMinutes);
    end
    pause(60*warmupTimeMinutes);
    if (verbose)
        fprintf('done.  Measuring.\n');
    end
    
    % Measure.
    for pp = 1:nPrimaries
        theScreenOnePrimarySettings = zeros(nPrimaries,1);
        theScreenOnePrimarySettings(pp) = 1;
        targetScreenSpdMeasured(:,pp) = MeasurePlainScreenSettings(theScreenOnePrimarySettings,...
                                     S,window,windowRect,'measurementOption',true,'verbose',verbose);
        clear theScreenOnePrimarySettings
        
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
% What we want is in targetScreenSpd, what we got is in
% targetScreenSpdMeasured.
% Plot the spd results.
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
if (MEASURE)
    screenCalObj = theData.screenCalObj;
    theData = rmfield(theData,'screenCalObj');
   
    %% Set the primaries in the calibration to the measured results.
    %
    % It's important to also set the sensor color space, because the
    % transform between sensor/primaries is cached when we set it.
    screenCalObj.set('P_device',targetScreenSpdMeasured);
    SetSensorColorSpace(screenCalObj,T_cones,S);
    
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
        % the background will be when we set the screen to the
        % background settings.  That then lets us compute contrast relative
        % to the background we're going to get. We want to do this from the
        % settings, so that the background isn't mucked up by quantization.
        screenBgSpd = PrimaryToSpd(screenCalObj,SettingsToPrimary(screenCalObj,theData.ptCldScreenSettingsCheckCal(:,1)));
        screenBgExcitations = T_cones * screenBgSpd;
        
        % Set up point cloud for finding best settings
        [contrastPtCld,ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',verbose);
        
        %% Generate some settings values corresponding to known contrasts
        %
        % The reason for this is to measure and check these.  This logic follows
        % how we handled an actual gabor image above. The quantization to
        % nQuantizeLevels isn't strictly needed, but nor is it doing harm.
        rawMonochromeUnquantizedContrastCheckCal = [0 0.25 -0.25 0.5 -0.5 1 -1];
        rawMonochromeContrastCheckCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
        desiredContrastCheckCal = theData.spatialGaborTargetContrast*theData.targetStimulusContrastDir*rawMonochromeContrastCheckCal;
        desiredExcitationsCheckCal = ContrastToExcitation(desiredContrastCheckCal,screenBgExcitations);
        
        % For each check calibration find the settings that
        % come as close as possible to producing the desired excitations.
        %
        % If we measure for a uniform field the spectra corresopnding to each of
        % the settings in the columns of ptCldScreenSettingsCheckCall, then
        % compute the cone contrasts with respect to the backgound (0 contrast
        % measurement, first settings), we should approximate the cone contrasts in
        % desiredContrastCheckCal.
        fprintf('Point cloud exhaustive method, finding settings\n')
        ptCldScreenSettingsCheckCal = SettingsFromPointCloud(contrastPtCld,desiredContrastCheckCal,ptCldSettingsCal);
        ptCldScreenPrimariesCheckCal = SettingsToPrimary(screenCalObj,ptCldScreenSettingsCheckCal);
        ptCldScreenSpdCheckCal = PrimaryToSpd(screenCalObj,ptCldScreenPrimariesCheckCal);
        ptCldScreenExcitationsCheckCal = SettingsToSensor(screenCalObj,ptCldScreenSettingsCheckCal);
        ptCldScreenContrastCheckCal = ExcitationsToContrast(ptCldScreenExcitationsCheckCal,screenBgExcitations);
    end
    
    %% Measure contrasts of the settings we computed in SpectralTestCal
    %
    % Measure the contrast points.  We've already got the settings so all we
    % need to do is loop through and set a uniform field to each of the
    % settings in ptCldScreenSettingsCheckCall and measure the corresponding
    % spd.
    [ptCldScreenSpdMeasuredCheckCal, ptCldScreenSettingsIntegersCheckCal] = MeasurePlainScreenSettings(ptCldScreenSettingsCheckCal,...
        S,window,windowRect,'measurementOption',true,'verbose',verbose);
   
end

%% Make plot of measured versus desired spds.
%
% The desired spds are in ptCldScreenSpdCheckCal
%
% Choose to use either raw or scaled measurement spectra.
whichToAnalyze = 'raw';
switch (whichToAnalyze)
    case 'raw'
        ptCldSpd = ptCldScreenSpdMeasuredCheckCal;
    case 'scaled'
        for tt = 1:nTestPoints
            % Find scale factor to best bring into alignment with target
            scaleSpdToTargetFactor(tt) = ptCldScreenSpdMeasuredCheckCal(:,tt)\ptCldScreenSpdCheckCal(:,pp);
            ptCldSpd(:,tt) = scaleSpdToTargetFactor(tt) * ptCldScreenSpdMeasuredCheckCal(:,tt);
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
    plot(wls,ptCldScreenSpdCheckCal(:,tt),'k-','LineWidth',3) % Target spectra
    plot(wls,ptCldSpd(:,tt),'r-','LineWidth',2); % Measured spectra
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
ptCldExcitationsMeasured = T_cones * ptCldSpd;
ptCldBgExcitationsMeasured = ptCldExcitationsMeasured(:,1);
ptCldScreenContrastMeasuredCheckCal = ExcitationToContrast(ptCldExcitationsMeasured,ptCldBgExcitationsMeasured);

% Add the nominal contrasts to be compared in the following graph. These
% should be all lined up on the 45-degree line in the following graph.
addNominalContrast = true;
if (addNominalContrast)
    ptCldExcitationsNominal = T_cones * ptCldScreenSpdCheckCal;
    ptCldBgExcitationsNominal = ptCldExcitationsNominal(:,1);
    ptCldContrastNominal = ExcitationToContrast(ptCldExcitationsNominal,ptCldBgExcitationsNominal);
else
end

% Plot measured versus desired contrasts
figure; hold on;
plot(desiredContrastCheckCal(1,:),ptCldScreenContrastMeasuredCheckCal(1,:),'ro','MarkerSize',14,'MarkerFaceColor','r');   % L - measured
plot(desiredContrastCheckCal(2,:),ptCldScreenContrastMeasuredCheckCal(2,:),'go','MarkerSize',12,'MarkerFaceColor','g');   % M - measured
plot(desiredContrastCheckCal(3,:),ptCldScreenContrastMeasuredCheckCal(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');   % S - measured
if (addNominalContrast)
    plot(desiredContrastCheckCal(1,:),ptCldContrastNominal(1,:),'ro','MarkerSize',19);   % L - target
    plot(desiredContrastCheckCal(2,:),ptCldContrastNominal(2,:),'go','MarkerSize',16);   % M - target
    plot(desiredContrastCheckCal(3,:),ptCldContrastNominal(3,:),'bo','MarkerSize',14);   % S - target
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

%% Close screen and save out the measurement data.
if (MEASURE)
    % Close
    CloseScreen;
    CloseSpectroradiometer;
    
    % Save data with the name containing dayTimestr, so that we can
    % automatically load in the most recent output.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,sprintf('testImageDataCheck_%s_%s',conditionName,dayTimestr));
        save(testFilename,'theData','targetScreenSpd','targetScreenSpdMeasured','screenCalObj', ...
            'screenBgSpd','screenBgExcitations','desiredContrastCheckCal','desiredExcitationsCheckCal', ...
            'ptCldScreenSettingsCheckCal','ptCldScreenPrimariesCheckCal','ptCldScreenSpdCheckCal','ptCldScreenExcitationsCheckCal','ptCldScreenContrastCheckCal', ...
            'ptCldScreenSettingsIntegersCheckCal','ptCldScreenSpdMeasuredCheckCal','ptCldScreenContrastMeasuredCheckCal', ...
            'ptCldExcitationsNominal','ptCldBgExcitationsNominal','ptCldContrastNominal');
    end
end
