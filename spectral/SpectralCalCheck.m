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

% In 'theData', following variables are available.
% ['S','T_cones','projectorCalObj','subprimaryCalObjs','projectorSettingsImage', ...
%        'projectorPrimaryPrimaries','projectorPrimarySettings','theDesiredContrastCheckCal', ...
%        'thePointCloudSettingsCheckCal','thePointCloudContrastCheckCal','thePointCloudSpdCheckCal']

% Set some variables.
S = theData.S; % Range of the spectrum.
wls = SToWls(S); % Wavelength. 
nPrimaries = 3; % Number of primaries.
nSubprimaries = theData.subprimaryCalObjs{1}.get('nDevices'); % Number of subprimaries.
subprimaryNInputLevels = size(theData.subprimaryCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:7 9:15];
nTestPoints = size(theData.thePointCloudContrastCheckCal,2);
T_cones = theData.T_cones;

%% Measure the desired target primaries
%
% This just needs to measure.
%
% Target Spds.
targetPrimarySpd = theData.projectorCalObj.get('P_device');

% Make a loop for measuring all primaries.
for pp = 1:nPrimaries
    [isolatingSpdMeasuredRaw(:,pp) isolatingSpdMeasuredNorm(:,pp), measuredToTarget(pp)] = MeasureDesiredTargetPrimaries(theData.projectorPrimaryPrimaries(:,pp), ...
                                               theData.subprimaryCalObjs{pp},pp,'projectorMode',true,'measurementOption',true,'verbose',false);
end

% Set the primaries as meausred results.
theData.projectorCalObj.set('P_device',isolatingSpdMeasuredRaw);

% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in targetPrimarySpd, what we got is in
% isolatingSpdMeasured.
% Plot the raw spd results.
figure; 
for pp = 1:nPrimaries
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetPrimarySpd(:,pp),'k-')
    plot(wls,isolatingSpdMeasuredRaw(:,pp),'r--');
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
end

% Plot the normalized spd results.
figure; 
for pp = 1:nPrimaries
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetPrimarySpd(:,pp),'k-')
    plot(wls,isolatingSpdMeasuredNorm(:,pp),'r--');
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
end

%% Set each primary to the settings we loaded in and measure
% Add VPixx toolbox ('Datapixx') to path.
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

%% Measure contrasts of the settings we computed in SpectralTestCal
%
% Measure the contrast points.  We've already got the settings so all we
% need to do is loop through and set a uniform field to each of the
% settings in thePointCloudSettingsCheckCal and measure the corresponding
% spd.
[thePointCloudSpdMeasured] = MeasureProjectorPrimarySettings(theData.thePointCloudSettingsCheckCal,...
    theData.projectorCalObj,T_cones,'projectorMode',true,'measurementOption',true,'verbose',true);

% thePointCloudSpdMeasured = MeasureSpdFromProjectorSettings(theData.thePointCloudSettingsCheckCal);

% Make plot of meeasured versus desired spds.  The desired spds are in
% theData.thePointCloudSpdCheckCal
for tt = 1:nTestPoints
    measuredToTarget(tt) = sum(theData.thePointCloudSpdCheckCal(:,tt))./sum(thePointCloudSpdMeasured(:,tt));
    thePointCloudSpdMeasuredNorm(:,tt) = thePointCloudSpdMeasured(:,tt) .* measuredToTarget(tt)
end

figure; 
for tt = 1:nTestPoints
    subplot(2,round(nTestPoints/2),tt); hold on; 
    plot(wls,theData.thePointCloudSpdCheckCal(:,tt),'k-') % Target spectra
    plot(wls,thePointCloudSpdMeasuredNorm(:,tt),'r--'); % Measured spectra
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured')
end

%% Compute cone contrasts for each spectrum relative to the background
%
% We use the fact that the background settings are in the first column of 
% theData.thePointCloudSettingsCheckCal.
testExcitations = T_cones * thePointCloudSpdMeasured;
bgExcitations = testExcitations(:,1);
testContrasts = (testExcitations - bgExcitations) ./ bgExcitations;

% Plot measured versus desired contrasts
figure; hold on;
plot(theData.thePointCloudContrastCheckCal(1,:),testContrasts(1,:),'r*'); % L
plot(theData.thePointCloudContrastCheckCal(2,:),testContrasts(2,:),'gO'); % M
plot(theData.thePointCloudContrastCheckCal(3,:),testContrasts(3,:),'b+'); % S
% axis('square');
xlabel('Desired contrast');
ylabel('Measured contrast');
% xlim([min(theData.thePointCloudContrastCheckCal,[],'all') max(theData.thePointCloudContrastCheckCal,[],'all')])
% ylim([min(theData.thePointCloudContrastCheckCal,[],'all') max(theData.thePointCloudContrastCheckCal,[],'all')])
legend('L','M','S','location','southeast')
axis('square')
title('Desired vs. Measured LMS Contrast')
