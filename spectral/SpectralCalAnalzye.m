% Read and analyze output of SpectralCalCheck
%

% History:
%    11/09/21  dhb  Wrote it.

%% Load output of SpectralCalCompute.
% 
% We can remove this now that we are saving testData from SpectralCalCheck,
% after next round of measurement.
conditionName = 'LminusMSmooth';
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theData = load(testFilename);
end

%% Load output of SpectralCalCheck
dayTimestr = '2021-11-08_15-43-59.mat';
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageDataCheck_%s_%s',conditionName,dayTimestr));
    theCheckData = load(testFilename);
    % theData = theCheckData.theData;
else
    error('Cannot find data file');
end

%% Set up some variables that we need
%
% This just needs to measure.
%
% Target Spds.
targetPrimarySpd = theData.projectorCalObj.get('P_device');

% Set some variables.
S = theData.S;                                                  % Range of the spectrum.
wls = SToWls(S);                                                % Wavelength. 
nPrimaries = size(targetPrimarySpd,2);                          % Number of primaries.
nSubprimaries = theData.subprimaryCalObjs{1}.get('nDevices');   % Number of subprimaries.
subprimaryNInputLevels = size(theData.subprimaryCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:7 9:15];
nTestPoints = size(theData.thePointCloudContrastCheckCal,2);
T_cones = theData.T_cones;

%% Look at measured primaries
%
% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in targetPrimarySpd, what we got is in
% isolatingSpdMeasured.
% Plot the spd results.
primarySpdFig = figure; clf; 
for pp = 1:nPrimaries
    checkTargetPrimaries(:,pp) = SettingsToPrimary(theData.subprimaryCalObjs{pp},theData.projectorPrimarySettings(:,pp));
    if (max(abs(checkTargetPrimaries(:,pp)-theData.projectorPrimaryPrimaries(:,pp))) > 0.005)
        error('Cannot reconstruct target primary values from target primary settings to quantization tolerance');
    end
    checkTargetPrimarySpd(:,pp) = PrimaryToSpd(theData.subprimaryCalObjs{pp},checkTargetPrimaries(:,pp));
    if (any(targetPrimarySpd(:,pp) ~= checkTargetPrimarySpd(:,pp)))
        error('Cannot reconstruct target primary spd from calibraion and primary values');
    end

    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetPrimarySpd(:,pp),'k','LineWidth',4)
    plot(wls,theCheckData.isolatingSpdMeasured(:,pp),'r','LineWidth',3);
    xlabel('Wavelength (nm)');
    ylabel('Spectral power distribution');
    title('Comparison of raw measured and desired spds');
end

%% Use regression to express measured primaries in terms of subprimary calibration
primaryValueFig = figure; clf; 
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize/3];
set(gcf,'position',figurePosition);
for pp = 1:nPrimaries
    P_device{pp} = theData.subprimaryCalObjs{pp}.get('P_device');
    P_ambient{pp} = theData.subprimaryCalObjs{pp}.get('P_ambient');
    if (any(P_ambient{pp}(:) ~= 0))
        error('Need to handle non-zero ambient');
    end
    regressPointCloudPrimaries(:,pp) = P_device{pp}\theCheckData.isolatingSpdMeasured(:,pp);

    figure(primaryValueFig);
    subplot(1,nPrimaries,pp); hold on;
    plot(checkTargetPrimaries(:,pp),regressPointCloudPrimaries(:,pp),'ro','MarkerSize',12,'MarkerFaceColor','r');
    xlim([0,1]); ylim([0,1]); axis('square');
    title('Regression versus desired primary values')
    xlabel('Desired Primary'); ylabel('Regression on Measured Primary');

    figure(primarySpdFig);
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,P_device{pp}*regressPointCloudPrimaries(:,pp),'g','LineWidth',2);
    legend('Target','Measured','Regress');
end

%% Make plot of measured versus desired spds.
%
% We should be able to obtain these as a linear combination
% of the primaries.
thePointCloudSpd = theCheckData.thePointCloudSpdMeasured;
figure;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize];
set(gcf,'position',figurePosition);
for tt = 1:nTestPoints
    % Find linear combination of measured primaries to produce test
    regressPointCloudPrimaries(:,tt) = theCheckData.isolatingSpdMeasured\theCheckData.thePointCloudSpdMeasured(:,tt);
    regressPointCloudSpd(:,tt) = theCheckData.isolatingSpdMeasured*regressPointCloudPrimaries(:,tt);
    subplot(round(nTestPoints/2),2,tt); hold on;
    plot(wls,theData.thePointCloudSpdCheckCal(:,tt),'k-','LineWidth',4) % Target spectra
    plot(wls,theCheckData.thePointCloudSpdMeasured(:,tt),'r-','LineWidth',3); % Measured spectra
    plot(wls,regressPointCloudSpd(:,tt),'g-','LineWidth',2); % Measured spectra
    xlabel('Wavelength (nm)')
    ylabel('Spectral power distribution')
    legend('Target','Measured','Regress')
    title(sprintf('Test %d raw',tt),'fontsize',16)
end

% Get primaries to compare regression primaries to
thePointCloudPrimariesCheckCal = SettingsToPrimary(theData.projectorCalObj,)
