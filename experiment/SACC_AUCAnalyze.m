% SACC_AUCAnalyze
%
% This is to make a summary plot of AUC of all subjects according to
% different filters.

% History:
%    5/2/23    smo    - Wrote it.
%    5/9/23    smo    - Added some more plots.

%% Initialize.
clear; close all;

%% Save the plots if you want.
SAVEPLOTS = true;
imgFileFormat = '.tiff';

%% Read the AUC summary table.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
    T = readtable(testFilename);
end

%% Calculate mean and standard error of AUC.
%
% Read and set some variables that we need.
filterOptions = unique(T.Filter);
subjectOptions = unique(T.Subject);
nFilters = length(filterOptions);
nSubjects = length(unique(T.Subject));
SFVals = [3 6 9 12 18];

% Extract AUC over the filters.
% The array of AUC is looks like Subjects (rows) x Filters (columns).
for ff = 1:nFilters
    AUC(:,ff) = T.AUC(strcmp(T.Filter, filterOptions(ff)));
end

% Calculate mean and standard error of AUC.
meanAUC = mean(AUC,1);
errorAUC = std(AUC,1)./sqrt(nSubjects);

%% Get log sensitivity per each filter.
for ff = 1:nFilters
    logSensitivity_3cpd(:,ff) = T.LogSensitivity_3cpd(strcmp(T.Filter, filterOptions(ff)));
    logSensitivity_6cpd(:,ff) = T.LogSensitivity_6cpd(strcmp(T.Filter, filterOptions(ff)));
    logSensitivity_9cpd(:,ff) = T.LogSensitivity_9cpd(strcmp(T.Filter, filterOptions(ff)));
    logSensitivity_12cpd(:,ff) = T.LogSensitivity_12cpd(strcmp(T.Filter, filterOptions(ff)));
    logSensitivity_18cpd(:,ff) = T.LogSensitivity_18cpd(strcmp(T.Filter, filterOptions(ff)));
end

% Log sensitivity. The array is sorted subject x spatial frequency x
% filter.
for ff = 1:nFilters
    logSensitivity(:,:,ff) = [logSensitivity_3cpd(:,ff) logSensitivity_6cpd(:,ff) logSensitivity_9cpd(:,ff) logSensitivity_12cpd(:,ff) logSensitivity_18cpd(:,ff)];
end

% Mean log sensitivity per each filter.
meanLogSensitivity = mean(logSensitivity);

%% Plotting the results from here.
%% 1) Plot AUC over the filters (bar graph)
figure; hold on;
xAxisTicks = [1:nFilters];
for ff = 1:nFilters
    bar(xAxisTicks(ff),meanAUC(ff),'facecolor',[ff*0.2 ff*0.2 0]);
end
xticks(xAxisTicks);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('Mean AUC','Fontsize',15);
errorbar(xAxisTicks,meanAUC,errorAUC,'k+','linewidth',1);
text(xAxisTicks, meanAUC/2, num2str(round(meanAUC,2)'),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'fontsize', 14, 'fontweight', 'bold', 'color', [0.75 0.75 0.75]);
title('Mean AUC results over the Filters','Fontsize',15);
subtitle(sprintf('Each bar is the average of (%d) subjects',nSubjects),'FontSize',13);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('AUC_Over_Filters',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 2) Plot AUC over the subjects (bar graph)
AUCSubFig = figure; clf; hold on;
figPosition = [100 100 1500 1000];
set(AUCSubFig, 'position', figPosition);
for ff = 1:nFilters
    subplot(2,3,ff);
    xAxisTicks = [1:nSubjects];
    bar(xAxisTicks,AUC(:,ff),'facecolor',[ff*0.2 ff*0.2 0]);
    xticks(xAxisTicks);
    xticklabels(subjectOptions);
    xlabel('Subject','Fontsize',15);
    ylabel('AUC','Fontsize',15);
    title(sprintf('Filter = %s',filterOptions{ff}),'fontsize',15);
    text(xAxisTicks(1),34,sprintf('Max = (%.2f) / Min = (%.2f) / Avg = (%.2f)', ...
        max(AUC(:,ff)), min(AUC(:,ff)), mean(AUC(:,ff))),...
        'fontsize',12);
end
sgtitle('AUC results over the subjects per each Filter','Fontsize',20);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('AUC_Over_Subjects',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 3) Plot AUC over the subjects and filters (line graph)
figure; clf; hold on;
xAxisTicks = [1:nFilters];
xticks(xAxisTicks);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('AUC','Fontsize',15);
for ss = 1:nSubjects
    plot(xAxisTicks,AUC(ss,:),'o-','linewidth',0.7);
end
plot(xAxisTicks,meanAUC,'k+-','color',[0 0 0 0.4],'linewidth',7);
text(xAxisTicks, meanAUC+0.3, num2str(round(meanAUC,2)'),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'fontsize', 11, 'fontweight', 'bold', 'color', 'k');
title('Mean AUC results over the Filters','Fontsize',15);
legendPlot = append('Sub ',subjectOptions);
legendPlot{end+1} = 'Mean';
legend(legendPlot,'location','northeastoutside');
title('AUC results over the filters and subjects','fontsize',15);
subtitle(sprintf('Mean result is the average of (%d) subjects',nSubjects),'FontSize',13);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('AUC_Over_SubjetsAndFilters',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 4-1) Plot Grand CCSF (line graph).
grandFig = figure; clf; hold on;
set(grandFig, 'position', figPosition);

for ff = 1:nFilters
    subplot(2,3,ff);
    plot(SFVals, meanLogSensitivity(:,:,ff), 'ko-','color',[ff*0.2 ff*0.2 0],'linewidth',3);
    xticks(SFVals);
    xticklabels(SFVals);
    xlabel('Spatial Frequency','Fontsize',15);
    ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', nSubjects),'Fontsize',15);
    yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
    ylim(log10([1 600]));
    yticks(yaxisRange);
    ytickformat('%.2f');
    legend(sprintf('Filter %s', filterOptions{ff}), 'fontsize', 15);
end
sgtitle(sprintf('Each figure is the averaged results of (%d) subjects', nSubjects), 'Fontsize',20);

% Grand grand mean.
meanLogSensitivityAllFilters = mean(squeeze(meanLogSensitivity),2);
subplot(2,3,6)
plot(SFVals, meanLogSensitivityAllFilters, 'k+-','color',[0 0 0 0.4],'linewidth',7);
xticks(SFVals);
xticklabels(SFVals);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', nSubjects*nFilters),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
legend('Mean of all filters', 'fontsize', 15);
title('This is the mean results of all filters and all subjects', 'fontsize',13);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('AUC_GrandMean',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 4-2) Plot Grand CCSF (line graph).
% We will draw the same graphs as above, but here we draw in one figure so
% that we can better see the difference over the filters.
figure; clf; hold on;

for ff = 1:nFilters
    plot(SFVals, meanLogSensitivity(:,:,ff), 'ko-','color',[ff*0.2 ff*0.2 0],'linewidth',2);
    xticks(SFVals);
    xticklabels(SFVals);
    xlabel('Spatial Frequency','Fontsize',15);
    ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', nSubjects),'Fontsize',15);
    yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
    ylim(log10([1 600]));
    yticks(yaxisRange);
    ytickformat('%.2f');    
end
legend(filterOptions, 'fontsize', 15);
title('Mean log sensitivity over spatial frequency per each filter','fontsize',15);


% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('AUC_GrandMean_Overlap',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end
