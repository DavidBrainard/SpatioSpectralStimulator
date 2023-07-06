% SACC_FitPFoverFocusLevel.
%
% This is to compare CS fitting results when doing experiment with and
% without setting a focus to the image by controlling the trombone.

% History:
%    07/05/23   smo   - Wrote it.

%% Initialize.
clear; close all;

%% Load data.
if (ispref('SpatioSpectralStimulator','SACCData'))
    subjectName = 'Semin';
    SF = '18_cpd';
    testFiledir = getpref('SpatioSpectralStimulator','SACCData');
    testFiledir = fullfile(testFiledir,subjectName,SF);
    fileList = dir(fullfile(testFiledir,'CS_Semin_*'));
end

% Read out the file names.
nFiles = length(fileList);
for ff = 1:nFiles
    filenameList{ff} = fileList(ff).name;
    theData(ff) = load(fullfile(testFiledir,filenameList{ff}));
end

%% Set variables.
PFMethod = 'weibull';
PF = @PAL_Weibull;
paramsFree = [1 1 0 1];
minSlope = 0.1;
maxSlope = 10;
nSlopes = 20;
slopeValList = 10.^linspace(log10(minSlope),log10(maxSlope),nSlopes);
axisLog = true;
VERBOSE = false;
bootConfInterval = 0.8;
focusLevels = {'focus','+0.5D','+1.5D','+2.5D'};

%% Fitting happens here.
for ff = 1:nFiles
    thresholdCriterion = 0.81606;
    [threshold, para, dataOut] = theData(ff).estimator.thresholdMLE('thresholdCriterion', thresholdCriterion, 'returnData', true);
    
    testContrastLog = dataOut.examinedContrasts;
    stimLevels = 10.^testContrastLog;
    pCorrect = dataOut.pCorrect;
    
    nTrials = theData(ff).estimator.nRepeat;
    
    % Fitting happens here.
    [paramsFitted,thresholdFitted] = FitPFToData(stimLevels, pCorrect, ...
        'PF', PFMethod, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
        'newFigureWindow', false, 'axisLog', axisLog,...
        'questPara', [],'addLegend',false, ...
        'beta',slopeValList,'nBootstraps',0,'bootConfInterval',bootConfInterval);
    
    % Get smooth fitted line.
    nFineStimLevels = 1000;
    fineStimLevels = linspace(0, max(stimLevels), nFineStimLevels);
    smoothPsychometric = PF(paramsFitted, fineStimLevels);
    
    % Switch to log space to plot it.
    stimLevelsPlot = log10(stimLevels);
    thresholdFittedLog = log10(thresholdFitted);
    fineStimLevelsPlot = log10(fineStimLevels);
    
    % Save out the data that needs to plot.
    stimLevelsPlotHandles{ff} = stimLevelsPlot;
    pCorrectHandles{ff} = pCorrect;
    fineStimLevelsPlotHandles{ff} = fineStimLevelsPlot;
    smoothPsychometricHandles{ff} = smoothPsychometric;
    thresholdFittedLogHandles{ff} = thresholdFittedLog;
end

%% Plot it.
%
% 1) One grand plot that contain everything.
figure;
figPosition = [0 0 1600 400];
set(gcf,'position',figPosition);

for ff = 1:nFiles
    subplot(1,4,ff); hold on;
    ylim([0 1]);
    xlim([-2.2 -1.2]);
    xlabel('Log contrast','fontsize',15);
    ylabel('pCorrect','fontsize',15);
    title(append('Focus level = ',focusLevels{ff}),'fontsize',15);
    subtitle(sprintf('Threshold = %.4f (linear: %.4f)', thresholdFittedLogHandles{ff}, 10^thresholdFittedLogHandles{ff}),'fontsize',15);
    
    % Raw data.
    h_data = scatter(stimLevelsPlotHandles{ff}, pCorrectHandles{ff}, 150,...
        'MarkerEdgeColor', zeros(1,3), 'MarkerFaceColor', ones(1,3) * 0.5, 'MarkerFaceAlpha', 0.5);
    
    % Smooth fit.
    h_pffit = plot(fineStimLevelsPlotHandles{ff},smoothPsychometricHandles{ff},'r','LineWidth',3);
    
    % Mark threshold.
    h_thresh = plot(thresholdFittedLogHandles{ff},thresholdCriterion,'ko','MarkerFaceColor','r','MarkerSize',12);
    legend('Data','PF-Fit','PF-Threshold','fontsize',13,'location','southwest');
end

% 2) Comparison of fitting graph.
figure; hold on;
ylim([0 1]);
xlim([-2.2 -1.2]);
xlabel('Log contrast','fontsize',15);
ylabel('pCorrect','fontsize',15);

lineTransparency = [1 0.6 0.4 0.15];
for ff = 1:nFiles
    plot(fineStimLevelsPlotHandles{ff},smoothPsychometricHandles{ff},'color',[1 0 0 lineTransparency(ff)], 'LineWidth',3);
    h_thresh = plot(thresholdFittedLogHandles{ff},thresholdCriterion,'ko','MarkerFaceColor','r','MarkerSize',12);
end

f = get(gca,'children');
legend(flip(f([2:2:8])),focusLevels,'fontsize',14,'location','southwest');
