% SpectralAnalyze
%
% Read and analyze output of SpectralCalCheck
%

% History:
%    11/09/21  dhb  Wrote it.

%% Initialize.
clear; close all;

%% Load output of SpectralCalCompute.
% 
% We can remove this now that we are saving testData from SpectralCalCheck,
% after next round of measurement.
conditionName = 'LminusMSmooth';
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = fullfile(testFiledir,sprintf('testImageData_%s',conditionName));
    theComputeData = load(testFilename);
end

%% Load output of SpectralCalCheck
if (ispref('SpatioSpectralStimulator','TestDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
    testFilename = GetMostRecentFileName(testFiledir,sprintf('testImageDataCheck_%s',conditionName),'olderDate',0);
    theCheckData = load(testFilename);
    theData = theCheckData.theData;
else
    error('Cannot find data file');
end

%% Set up some variables that we need
%
% Target Spds.  The check here just verifies that we are consistent between
% current compute code and this measurement.
targetScreenSpd = theCheckData.targetScreenSpd;
targetScreenSpdCompute = theComputeData.screenCalObj.get('P_device');
if (any(targetScreenSpd(:) ~= targetScreenSpdCompute(:)))
    error('Strange change in target spd');
end
%targetScreenSpd = theComputeData.screenCalObj.get('P_device');

% Set some variables.
S = theData.S;                                                  % Range of the spectrum.
wls = SToWls(S);                                                % Wavelength. 
nPrimaries = size(targetScreenSpd,2);                          % Number of primaries.
nSubprimaries = theData.subprimaryCalObjs{1}.get('nDevices');   % Number of subprimaries.
subprimaryNInputLevels = size(theData.subprimaryCalObjs{1}.get('gammaInput'),1);
logicalToPhysical = [0:7 9:15];
nTestPoints = size(theData.thePtCldScreenContrastCheckCal,2);
T_cones = theData.T_cones;

%% Look at measured primaries
%
% Make plot comparing what we wanted for primaries versus what we got.
% What we want is in targetScreenSpd, what we got is in
% targetScreenSpdMeasured.
% Plot the spd results.
primarySpdFig = figure; clf; 
for pp = 1:nPrimaries
    checkTargetPrimaries(:,pp) = SettingsToPrimary(theData.subprimaryCalObjs{pp},theData.projectorPrimarySettings(:,pp));
    if (max(abs(checkTargetPrimaries(:,pp)-theData.projectorPrimaryPrimaries(:,pp))) > 0.005)
        error('Cannot reconstruct target primary values from target primary settings to quantization tolerance');
    end
    checkTargetSettings(:,pp) = PrimaryToSettings(theData.subprimaryCalObjs{pp},theData.projectorPrimaryPrimaries(:,pp));
    if (any(theData.projectorPrimarySettings(:,pp) ~= checkTargetSettings(:,pp)))
        error('Cannot reconstruct target primary settings from calibraion and primary values');
    end
    checktargetScreenSpd(:,pp) = PrimaryToSpd(theData.subprimaryCalObjs{pp},checkTargetPrimaries(:,pp));
    if (any(targetScreenSpd(:,pp) ~= checktargetScreenSpd(:,pp)))
        error('Cannot reconstruct target primary spd from calibraion and primary values');
    end

    subplot(nPrimaries,1,pp); hold on;
    plot(wls,targetScreenSpd(:,pp),'k','LineWidth',4)
    plot(wls,theCheckData.targetScreenSpdMeasured(:,pp),'r','LineWidth',3);
    plot(wls,checktargetScreenSpd(:,pp),'g','LineWidth',1);
    legend('Target','Measured','Target Check');
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
    regressTargetPrimaries(:,pp) = P_device{pp}\theCheckData.targetScreenSpdMeasured(:,pp);

    figure(primaryValueFig);
    subplot(1,nPrimaries,pp); hold on;
    plot(checkTargetPrimaries(:,pp),regressTargetPrimaries(:,pp),'ro','MarkerSize',12,'MarkerFaceColor','r');
    xlim([0,1]); ylim([0,1]); axis('square');
    plot([0 1],[0 1],'k');
    title('Regression versus desired subprimary values')
    xlabel('Desired Primary'); ylabel('Regression on Measured Primary');

    figure(primarySpdFig);
    subplot(nPrimaries,1,pp); hold on;
    plot(wls,P_device{pp}*regressTargetPrimaries(:,pp),'b','LineWidth',1);
    legend('Target','Measured','Target Check','Regress');
end

%% Make plot of measured versus desired spds.
%
% We should be able to obtain these as a linear combination
% of the primaries.
plotRegress = false;
figure; clf;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize];
set(gcf,'position',figurePosition);
for tt = 1:nTestPoints
    % Find linear combination of measured primaries to produce test
    regressPointCloudPrimaries(:,tt) = theCheckData.targetScreenSpdMeasured\theCheckData.thePtCldScreenSpdMeasuredCheckCal(:,tt);
    regressPointCloudSpd(:,tt) = theCheckData.targetScreenSpdMeasured*regressPointCloudPrimaries(:,tt);
    
    subplot(round(nTestPoints/2),2,tt); hold on;
    plot(wls,theCheckData.thePtCldScreenSpdCheckCal(:,tt),'k-','LineWidth',4)  % Target spectra
    plot(wls,theCheckData.thePtCldScreenSpdMeasuredCheckCal(:,tt),'r-','LineWidth',2); % Measured spectra
    xlabel('Wavelength (nm)');
    ylabel('Spectral power distribution');
    title(sprintf('Test %d raw',tt),'fontsize',16);
    if (plotRegress)
        plot(wls,regressPointCloudSpd(:,tt),'g-','LineWidth',1);              % Regression fit
        legend('Target','Measured','Regress');
    else
        legend('Target','Measured');
    end
end

% Another way of comparing measured and nominal spectra
thePointCloudSpdScatter = theCheckData.thePtCldScreenSpdMeasuredCheckCal;
figure; clf;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize];
set(gcf,'position',figurePosition);
spdLim = 2.5e-3;
for tt = 1:nTestPoints
    scaleFactor(tt) = theCheckData.thePtCldScreenSpdCheckCal(:,tt)\theCheckData.thePtCldScreenSpdMeasuredCheckCal(:,tt);
    subplot(round(nTestPoints/2),2,tt); hold on;
    plot(theCheckData.thePtCldScreenSpdCheckCal(:,tt),theCheckData.thePtCldScreenSpdMeasuredCheckCal(:,tt),'r+');
    xlabel('Nominal Spd');
    ylabel('Measured Spd');
    xlim([0 spdLim]); ylim([0 spdLim]);
    axis('square');
    title(sprintf('Test %d raw, factor %0.3f',tt,scaleFactor(tt)),'fontsize',12)
end


% Here we will compare the primaries we get from regression above to those
% we wanted, and see how we do.
figure; clf;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize/3];
set(gcf,'position',figurePosition);
for pp = 1:nPrimaries
    subplot(1,nPrimaries,pp); hold on;
    plot(theCheckData.thePtCldScreenPrimariesCheckCal(pp,:),regressPointCloudPrimaries(pp,:),'ro','MarkerSize',12,'MarkerFaceColor','r');
    plot([0 1],[0 1],'k');
    xlim([0,1]); ylim([0,1]); axis('square');
    title('Regression versus desired primary values')
    xlabel('Desired Projector Primary'); ylabel('Regression on Measured Projector Primary');
end

%% Explicitly compute some contrasts
analyzePointCloudExcitationsMeasured = T_cones * theCheckData.thePtCldScreenSpdMeasuredCheckCal;
analyzePointCloudBgExcitationsMeasured = analyzePointCloudExcitationsMeasured(:,1);
analyzePointCloudContrastMeasured = ExcitationToContrast(analyzePointCloudExcitationsMeasured,analyzePointCloudBgExcitationsMeasured);
if (any(analyzePointCloudContrastMeasured(:) ~= theCheckData.thePtCldScreenContrastMeasuredCheckCal(:)))
    error('Cannot get same measured contrasts in two different places');
end

analyzePointCloudExcitationsNominal = T_cones * theCheckData.thePtCldScreenSpdCheckCal;
analyzePointCloudBgExcitationsNominal = analyzePointCloudExcitationsNominal(:,1);
analyzePointCloudContrastNominal = ExcitationToContrast(analyzePointCloudExcitationsNominal,analyzePointCloudBgExcitationsNominal);
if (any(analyzePointCloudContrastNominal(:) ~= theCheckData.thePointCloudContrastNominal(:)))
    error('Cannot get same nominalcontrasts in two different places');
end

%% Plot measured versus desired contrasts
contrastFig = figure; hold on;
figureSize = 1000;
figurePosition = [1200 300 figureSize figureSize/3];
set(gcf,'position',figurePosition);

axisLim = 0.05;
theColors = ['r' 'g' 'b'];
for pp = 1:nPrimaries
    subplot(1,nPrimaries,pp); hold on;
    plot(theCheckData.theDesiredContrastCheckCal(pp,:),theCheckData.thePtCldScreenContrastMeasuredCheckCal(pp,:),[theColors(pp) 'o'],'MarkerSize',14,'MarkerFaceColor',theColors(pp));
    plot(theCheckData.theDesiredContrastCheckCal(pp,:),theCheckData.thePointCloudContrastNominal(pp,:), [theColors(pp) 'o'],'MarkerSize',18);
    plot(theCheckData.theDesiredContrastCheckCal(pp,1),theCheckData.thePtCldScreenContrastMeasuredCheckCal(pp,1),'ko','MarkerSize',14,'MarkerFaceColor','k');
    plot(theCheckData.theDesiredContrastCheckCal(pp,1),theCheckData.thePointCloudContrastNominal(pp,1), 'ko','MarkerSize',18);

    plot([-1 1],[-1 1],'k');
    xlabel('Desired contrast');
    ylabel('Measured contrast');
    xlim([-axisLim axisLim]);
    ylim([-axisLim axisLim]);
    axis('square');
    xlabel('Desired contrast');
    ylabel('Measured contrast');
    legend({'Measured','Nominal'},'location','southeast');
    title(sprintf('Cone class %d',pp));
end
