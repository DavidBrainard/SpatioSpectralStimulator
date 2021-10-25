% SpectralCalCheck
%
% Read in output of SpectralTestCal and make some measurements.
% 
% 10/19/2021  dhb,smo  Started on it

%% Initialize.
clear; close all;

%% Load output of SpectralTestCal
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,'testImageData1');
    theData = load(testFilename);
end

%% Measure the desired target primaries
%
% This just needs to measure.
%
% Target Spds.
targetPrimarySpd = theData.projectorCalObj.get('P_device');

% In 'theData', following variables are available.
% ['S','T_cones','projectorCalObj','subprimaryCalObjs','projectorSettingsImage', ...
%        'projectorPrimaryPrimaries','projectorPrimarySettings','theDesiredContrastCheckCal', ...
%        'thePointCloudSettingsCheckCal','thePointCloudContrastCheckCal','thePointCloudSpdCheckCal']

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
MEASURE = false;
if (MEASURE)
    for pp = 1:nPrimaries
        isolatingSpdMeasured(:,pp) = MeasureDesiredTargetPrimaries(theData.projectorPrimaryPrimaries(:,pp), ...
            theData.subprimaryCalObjs{pp},pp,'projectorMode',true,'measurementOption',true,'verbose',true);
    end
else
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = fullfile(testFiledir,'testImageData1Check');
        load(testFilename,'isolatingSpdMeasured','thePointCloudSpdMeasured','testContrasts');
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
    plot(wls,targetPrimarySpd(:,pp),'k-')
    plot(wls,isolatingSpdMeasured(:,pp),'r--');
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
    plot(wls,targetPrimarySpd(:,pp),'k-')
    plot(wls,scalePrimaryToTargetFactor(pp)*isolatingSpdMeasured(:,pp),'r--');
    xlabel('Wavelength (nm)');
    ylabel('Spectral power distribution');
    legend('Target','Measured');
    title('Comparison of raw measured and desired spds');
end

%% Set each primary to the settings we loaded in and measure - NEED TO MODIFY TO USE TbTb
%
% Add VPixx toolbox ('Datapixx') to path.
if (MEASURE)
    addpath(genpath('/home/colorlab/Documents/MATLAB/toolboxes/VPixx'));
    
    % Connect to the projector.
    isReady = Datapixx('open');
    isReady = Datapixx('IsReady');
    
    % Set the projector subprimaries here.
    for ss = 1:nSubprimaries
        Datapixx('SetPropixxHSLedCurrent', 0, logicalToPhysical(ss), round(theData.projectorPrimarySettings(ss,1)*(subprimaryNInputLevels-1))); % Primary 1
        Datapixx('SetPropixxHSLedCurrent', 1, logicalToPhysical(ss), round(theData.projectorPrimarySettings(ss,2)*(subprimaryNInputLevels-1))); % Primary 2
        Datapixx('SetPropixxHSLedCurrent', 2, logicalToPhysical(ss), round(theData.projectorPrimarySettings(ss,3)*(subprimaryNInputLevels-1))); % Primary 3
    end
    
    %% Set the primaries in the calibration to the meausred results.
    theData.projectorCalObj.set('P_device',meanPrimaryScaleFactor*isolatingSpdMeasured);
    
    %% Optional recompute of target settings
    %
    % We just measured the primaries.  So, we could recompute the contrast
    % check settings based on these measured primaries, rather than the
    % nominal ones we used in SpectralTestCal.  Here is where that happens
    RECOMPUTE = true;
    if (RECOMPUTE)
        % Parameters.  Should save these in SpectralTestCal and use here.
        nQuantizeBits = 14;
        nQuantizeLevels = 2^nQuantizeBits;
        projectorNInputLevels = 256;
        targetStimulusContrastDir = [1 -1 0]'; targetStimulusContrastDir = targetStimulusContrastDir/norm(targetStimulusContrastDir);
        spatialGaborTargetContrast = 0.04;

        % Figure out desired background excitations
        projectorBgExcitations = T_cones * theData.thePointCloudSpdCheckCal(:,1);
        
        % Build point cloud
        tic;
        fprintf('Point cloud exhaustive method, setting up cone contrast cloud, this takes a while\n')
        allProjectorSettingsCal = zeros(3,projectorNInputLevels^3);
        idx = 1;
        for ii = 0:(projectorNInputLevels-1)
            for jj = 0:(projectorNInputLevels-1)
                for kk = 0:(projectorNInputLevels-1)
                    allProjectorSettingsCal(:,idx) = [ii jj kk]'/projectorNInputLevels;
                    idx = idx+1;
                end
            end
        end

        % Get LMS excitations for each triplet of projector settings, and build a
        % point cloud object from these.
        allProjectorExcitations = SettingsToSensor(theData.projectorCalObj,allProjectorSettingsCal);
        allProjectorContrast = ExcitationsToContrast(allProjectorExcitations,projectorBgExcitations);
        allSensorPtCloud = pointCloud(allProjectorContrast');
        
        % Force point cloud setup by finding one nearest neighbor. This is slow,
        % but once it is done subsequent calls are considerably faster.
        findNearestNeighbors(allSensorPtCloud,[0 0 0],1);
        toc
        
        %% Generate some settings values corresponding to known contrasts % (THIS PART MAY BE GOING TO BE IN A FUNCTION LATER ON - SEMIN)
        %
        % The reason for this is to measure and check these.  This logic follows
        % how we handled an actual gabor image above.
        rawMonochromeUnquantizedContrastCheckCal = [0 0.25 -0.25 0.5 -0.5 1 -1];
        rawMonochromeContrastGaborCal = 2*(PrimariesToIntegerPrimaries((rawMonochromeUnquantizedContrastCheckCal +1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
        theDesiredContrastCheckCal = spatialGaborTargetContrast*targetStimulusContrastDir*rawMonochromeContrastGaborCal;
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
        theData.thePointCloudSettingsCheckCal = zeros(3,size(theDesiredContrastCheckCal,2));
        for ll = 1:size(theDesiredContrastCheckCal,2)
            minIndex = findNearestNeighbors(allSensorPtCloud,theDesiredContrastCheckCal(:,ll)',1);
            theData.thePointCloudSettingsCheckCal(:,ll) = allProjectorSettingsCal(:,minIndex);
        end
        theData.thePointCloudPrimariesCheckCal = SettingsToPrimary(theData.projectorCalObj,theData.thePointCloudSettingsCheckCal);
        theData.thePointCloudSpdCheckCal = PrimaryToSpd(theData.projectorCalObj,theData.thePointCloudPrimariesCheckCal);
        theData.thePointCloudExcitationsCheckCal = SettingsToSensor(theData.projectorCalObj,theData.thePointCloudSettingsCheckCal);
        theData.thePointCloudContrastCheckCal = ExcitationsToContrast(theData.thePointCloudExcitationsCheckCal,projectorBgExcitations);
      
    end
    
    %% Measure contrasts of the settings we computed in SpectralTestCal - SEE IF WE CAN WRITE A MEAUSURE SETTINGS ROUTINE THAT DOESN'T NEED CAL FILE OR T_CONES
    %
    % Measure the contrast points.  We've already got the settings so all we
    % need to do is loop through and set a uniform field to each of the
    % settings in thePointCloudSettingsCheckCal and measure the corresponding
    % spd.
    [thePointCloudSpdMeasured] = MeasureProjectorPrimarySettings(theData.thePointCloudSettingsCheckCal,...
        theData.projectorCalObj,T_cones,'projectorMode',true,'measurementOption',true,'verbose',true);
    % thePointCloudSpdMeasured = MeasureSpdFromProjectorSettings(theData.thePointCloudSettingsCheckCal);
end

%% Make plot of meeasured versus desired spds.
%
% The desired spds are in theData.thePointCloudSpdCheckCal

% Raw spds
figure; 
for tt = 1:nTestPoints
    subplot(2,round(nTestPoints/2),tt); hold on; 
    plot(wls,theData.thePointCloudSpdCheckCal(:,tt),'k-') % Target spectra
    plot(wls,thePointCloudSpdMeasured(:,tt),'r--'); % Measured spectra
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
    title('Raw Spds')
end

% Scale measured spds to target and then plot
for tt = 1:nTestPoints
    % Find scale factor to best bring into alignment with target
    scaleSpdToTargetFactor(tt) = thePointCloudSpdMeasured(:,tt)\theData.thePointCloudSpdCheckCal(:,pp);
    thePointCloudSpdMeasuredScaled(:,tt) = scaleSpdToTargetFactor(tt)*thePointCloudSpdMeasured(:,tt);
end
figure; 
for tt = 1:nTestPoints
    subplot(2,round(nTestPoints/2),tt); hold on; 
    plot(wls,theData.thePointCloudSpdCheckCal(:,tt),'k-') % Target spectra
    plot(wls,thePointCloudSpdMeasuredScaled(:,tt),'r--'); % Measured spectra
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
    title('Normalized Spds')
end

%% Compute cone contrasts for each spectrum relative to the background
%
% We use the fact that the background settings are in the first column of 
% theData.thePointCloudSettingsCheckCal.
whichToAnalyze = 'scaled';
switch (whichToAnalyze)
    case 'raw'
        testExcitations = T_cones * thePointCloudSpdMeasured;
    case 'scaled'
        testExcitations = T_cones * thePointCloudSpdMeasuredScaled;
    otherwise
        error('Unknown analyze type specified');
end
bgExcitations = testExcitations(:,1);
testContrasts = (testExcitations - bgExcitations) ./ bgExcitations;

% Add the target contrasts to be compared in the following graph.
addTargetContrast = true;
if (addTargetContrast)
    targetExcitations = T_cones * theData.thePointCloudSpdCheckCal;
    targetBgExcitations = targetExcitations(:,1);
    targetContrasts = (targetExcitations - targetBgExcitations) ./ targetBgExcitations;
else
    targetContrasts = zeros(1,size(testContrasts,2)); % If not passing this part, make it negative value not to be seen on the following graph.
end

% Plot measured versus desired contrasts
figure; hold on;
% Measured contrasts.
plot(theData.thePointCloudContrastCheckCal(1,:),testContrasts(1,:),'ro','MarkerSize',14,'MarkerFaceColor','r');   % L - measured
plot(theData.thePointCloudContrastCheckCal(2,:),testContrasts(2,:),'go','MarkerSize',12,'MarkerFaceColor','g');   % M - measured
plot(theData.thePointCloudContrastCheckCal(3,:),testContrasts(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');   % S - measured
% Target contrasts.
if (addTargetContrast)
    plot(theData.thePointCloudContrastCheckCal(1,:),targetContrasts(1,:),'ro','MarkerSize',19);   % L - target
    plot(theData.thePointCloudContrastCheckCal(2,:),targetContrasts(2,:),'go','MarkerSize',16);   % M - target
    plot(theData.thePointCloudContrastCheckCal(3,:),targetContrasts(3,:),'bo','MarkerSize',14);   % S - target
end
xlabel('Desired contrast');
ylabel('Measured contrast');
axisLim = 0.05;
xlim([-axisLim axisLim]);
ylim([-axisLim axisLim]);
axis('square');
line([-axisLim,axisLim], [-axisLim,axisLim], 'LineWidth', 1, 'Color', 'k');
grid on;
if (addTargetContrast)
    legend('L-Test','M-Test','S-Test','L-Target','M-Target','S-Target','location','southeast');
else
    legend('L','M','S','location','southeast');
end
title(sprintf('Desired vs. Measured LMS Contrast, %s',whichToAnalyze));

%% Save out the measurement data.
if (MEASURE)
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = fullfile(testFiledir,'testImageData1Check');
        save(testFilename,'isolatingSpdMeasured','thePointCloudSpdMeasured','testContrasts');
    end
end