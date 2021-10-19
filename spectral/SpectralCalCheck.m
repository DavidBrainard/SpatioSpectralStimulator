% SpectralCalCheck
%
% Read in output of SpectralTestCal and make some measurements.

%% Load output of SpectralTestCal
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,'testImageData1');
    theData = load(testFilename);
end
S = theData.S;

%% Set each primary to the settings we loaded in and measure
%
% Set the primaries here.  The primary settings are in the columns of
% theData.projectorPrimarySettings.

%% Measure the desired target primaries
%
% This just needs to measure.
%
% Make a loop for measuring all primaries.
targetPrimarySpd = projectorCalObj.get('P_device');
subprimaryNInputLevels = size(heData.projectorPrimaryPrimaries,1);
for pp = 1:nPrimaries
    isolatingSpdMeasured(:,pp) = MeasureDesiredTargetPrimaries(theData.projectorPrimaryPrimaries(:,pp),subprimaryNInputLevels,theData.subprimaryCalObjs{pp},pp,'projectorMode',true,'measurementOption',false,'verbose',false);
end
projectorCalObj.set('P_device',isolatingSpdMeasured);

% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in targetPrimarySpd, what we got is in
% isolatingSpdMeasured.

%% Measure contrasts of the settings we computed in SpectralTestCal
%
% Measure the contrast points.  We've already got the settings so all we
% need to do is loop through and set a uniform field to each of the
% settings in thePointCloudSettingsCheckCal and measure the corresponding
% spd.
% [thePointCloudSpdMeasured,projectorBgSpdMeasured] = MeasureLMSContrastGaborPatch_copy(thePointCloudSettingsCheckCal,projectorBgSettings,...
%     projectorCalObj,subprimaryCalObjs,T_cones,subprimaryNInputLevels,'projectorMode',true,'measurementOption',true,'verbose',true);
thePointCloudSpdMeasured = MeasureSpdFromProjectorSettings(theData.thePointCloudSettingsCheckCal);

% Make plot of meeasured versus desired spds.  The desired spds are in
% theData.thePointCloudSpdCheckCal

%% Compute cone contrasts for each spectrum relative to the background
%
% We use the fact that the background settings are in the first column of 
% theData.thePointCloudSettingsCheckCal.
testExcitations = theData.T_cones*measure.testSpd;
bgExcitations = testExcitations(:,1);
testContrasts = (testExcitations-bgExcitations) ./ bgExcitations;

% Plot measured versus desired contrasts
figure; hold on;
plot(thePointCloudContrastCheckCal(1,:),testContrasts(1,:),'r+'); % L
plot(thePointCloudContrastCheckCal(2,:),testContrasts(2,:),'g+'); % M
plot(thePointCloudContrastCheckCal(3,:),testContrasts(3,:),'b+'); % S
axis('square');
xlabel('Desired contrast');
ylabel('Measured contrast');