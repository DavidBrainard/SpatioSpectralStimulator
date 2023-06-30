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
SFOptionsStr = {'3 cpd', '6 cpd', '9 cpd', '12 cpd', '18 cpd'};
SFOptions = [3 6 9 12 18];
nSFs = length(SFOptions);

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
for ss = 1:nSFs
    logSensitivityAllFiltersTemp = squeeze(logSensitivity(:,ss,:));
    logSensitivityAllFiltersTemp = reshape(logSensitivityAllFiltersTemp, [size(logSensitivityAllFiltersTemp,1)*size(logSensitivityAllFiltersTemp,2) 1]);
    stdErrorAllFiltersAndSubjects(ss) = std(logSensitivityAllFiltersTemp)/sqrt(length(logSensitivityAllFiltersTemp));
end

%% Get subject gender info, MPOD, Iso-luminance results.
%
% Subject gender.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'Experiment_Progress.xlsx');
    tableGender = readtable(testFilename,'UseExcel',true,'sheet','SACC_Exp_Progress','range','B2:C48');
end

% MPOD.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'MPOD.mat');
    dataMPOD = load(testFilename);
end
MPOD = dataMPOD.MPOD;

% Iso-luminance determination results.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'IsoLuminanceDetermination.mat');
    dataFlicker = load(testFilename);
end
Flicker = dataFlicker.Flicker;

%% Plotting the results from here.
%% 1) Plot AUC over the filters (bar graph)
figure; hold on;
xAxisTicks = [1:nFilters];
for ff = 1:nFilters
    bar(xAxisTicks(ff),meanAUC(ff),'facecolor',[0.2 0.2 0.2]);
end
xticks(xAxisTicks);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('Mean AUC','Fontsize',15);
errorbar(xAxisTicks,meanAUC,errorAUC,'k+','linewidth',1);
text(xAxisTicks, meanAUC/2, num2str(round(meanAUC,2)'),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'fontsize', 14, 'fontweight', 'bold', 'color', [1 1 1]);
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
    plot(SFOptions, meanLogSensitivity(:,:,ff), 'ko-','color',[ff*0.2 ff*0.2 0],'linewidth',3);
    xticks(SFOptions);
    xticklabels(SFOptions);
    xlabel('Spatial Frequency','Fontsize',15);
    ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', nSubjects),'Fontsize',15);
    yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
    ylim(log10([1 600]));
    yticks(yaxisRange);
    ytickformat('%.2f');
    legend(sprintf('Filter %s', filterOptions{ff}), 'fontsize', 15);
end
sgtitle(sprintf('Each figure is the averaged results of (%d) subjects', nSubjects), 'Fontsize',20);

%% Grand grand mean.
meanLogSensitivityAllFilters = mean(squeeze(meanLogSensitivity),2);
% subplot(2,3,6); hold on;
figure; hold on;

% Plot error bar and delete the lines between adjacent points.
errorbar(SFOptions, meanLogSensitivityAllFilters, std(squeeze(meanLogSensitivity)')./sqrt(nFilters), 'k-');
set(get(gca, 'Children'), 'LineStyle', 'none');

% Plot the mean sensitivity data.
plot(SFOptions, meanLogSensitivityAllFilters, 'k.-', 'markersize',5, 'color',[0 0 0 0.4],'linewidth',4);
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', nSubjects*nFilters),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
xlim([min(SFOptions)-2 max(SFOptions)+2]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
legend('SEM of the filters (N=5)', 'Mean of all filters and subjects (N=160)', 'fontsize', 12);
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
    plot(SFOptions, meanLogSensitivity(:,:,ff), 'ko-','color',[ff*0.2 ff*0.2 0],'linewidth',2);
    xticks(SFOptions);
    xticklabels(SFOptions);
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

%% 5-1) Sensitivity difference across the filters.
% As a function of log sensitivity (filter A).
logSensitivityFilterA = logSensitivity(:,:,1);

% Loop over the filters.
for ff = 2:nFilters
    logSensitivityFilterTest = logSensitivity(:,:,ff);
    
    % Temp option for Semin's Laptop.
    figPosition = [0 0 1000 1000];
    fontSize = 10;
    
    % Figure to plot it as a function of log sensitivity (filter A).
    fig = figure; clf; hold on;
    set(fig,'position', figPosition);
    
    % Loop over spatial frequency.
    for ss = 1:nSFs
        subplot(2,3,ss); hold on;
        
        % Calculate error.
        logSensitivityError = logSensitivityFilterTest(:,ss)-logSensitivityFilterA(:,ss);
        meanLogSensitivityError = mean(logSensitivityError);
        stdErrorLogSensitivityError = std(logSensitivityError)/sqrt(length(logSensitivityError));
        
        % Sensitivity difference data.
        plot(logSensitivityFilterA(:,ss), logSensitivityError, 'o',...
            'markeredgecolor','k','markersize',6,'markerfacecolor',[ff*0.2 ff*0.2 0]);
        
        % Mean sensitivity line.
        plot([-3 +3], [meanLogSensitivityError meanLogSensitivityError], 'b-','linewidth',2,'color',[0 0 0.8 0.7]);
        
        % Standard error line.
        stdError_pos = meanLogSensitivityError + stdErrorLogSensitivityError;
        stdError_neg = meanLogSensitivityError - stdErrorLogSensitivityError;
        plot([-3 +3], [stdError_pos stdError_pos], 'k--','linewidth',1,'color',[0 0 0.8 0.9]);
        plot([-3 +3], [stdError_neg stdError_neg], 'k--','linewidth',1,'color',[0 0 0.8 0.9]);
        
        % No difference line for the reference.
        plot([-3 +3], [0 0], 'k-','linewidth',5,'color',[0.3 0.3 0.3 0.4]);
        
        xlabel('Log sensitivity (Filter A)','fontsize',fontSize);
        ylabel(sprintf('Log sensitivity Difference (Filter %s-A)',filterOptions{ff}),'fontsize',fontSize);
        xlim([1 2.6]);
        ylim([-0.6 0.6]);
        yticks([-0.6:0.2:0.6]);
        yticklabels([-0.6:0.2:0.6]);
        title(SFOptionsStr{ss},'fontsize',fontSize);
        text(1.05,0.55,sprintf('Max abs diffrence = %.2f (log) / %.2f (linear)', max(abs(logSensitivityError)), 10^max(abs(logSensitivityError))),'fontsize',fontSize-2);
        text(1.01,0.04,'No difference line','fontsize',fontSize-2);
        legend('Individual data', sprintf('Mean (N=%d)', nSubjects),'location','southeast','fontsize',fontSize-3);
    end
    sgtitle(sprintf('Log sensitivity difference between Filter ( A ) and Filter ( %s )',filterOptions{ff}),'Fontsize',fontSize+5);
    
    % Save the plot.
    if (SAVEPLOTS)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
            testFilename = fullfile(testFiledir, append(sprintf('SensitivityDifference_A_and_%s',filterOptions{ff}),imgFileFormat));
            saveas(gcf,testFilename);
            fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
        end
    end
end

%% 5-2) Sensitivity difference across the filters.
% As a function of spatial frequency.
logSensitivityFilterA = logSensitivity(:,:,1);

% Loop over the filters.
for ff = 2:nFilters
    logSensitivityFilterTest = logSensitivity(:,:,ff);
    
    % Temp option for Semin's Laptop.
    fontSize = 10;
    
    % Figure to plot it as a function of log sensitivity (filter A).
    figure; clf; hold on;
    
    % Loop over spatial frequency.
    for ss = 1:nSubjects
        % Calculate error.
        logSensitivityError = logSensitivityFilterTest(ss,:)-logSensitivityFilterA(ss,:);
        meanLogSensitivityError = mean(logSensitivityError);
        stdErrorLogSensitivityError = std(logSensitivityError)/sqrt(length(logSensitivityError));
        
        % Sensitivity difference data.
        plot(SFOptions, logSensitivityError, 'o-','markersize',6);
        
        xlabel('Spatial Frequency (cpd)','fontsize',fontSize);
        ylabel(sprintf('Log sensitivity Difference (Filter %s-A)',filterOptions{ff}),'fontsize',fontSize);
        xlim([min(SFOptions)-1 max(SFOptions)+1]);
        xticks(SFOptions);
        xticklabels(SFOptions);
        ylim([-0.6 0.6]);
        yticks([-0.6:0.2:0.6]);
        yticklabels([-0.6:0.2:0.6]);
    end
    
    % Mean results.
    plot(SFOptions, mean(logSensitivityFilterTest-logSensitivityFilterA), 'k+-', 'color', [0 0 1 0.5], ...
        'markerfacecolor', 'b', 'markeredgecolor', 'b', 'markersize', 5, 'linewidth', 5);
    
    % No difference line.
    plot([min(SFOptions)-1 max(SFOptions+1)], [0 0], 'k-','linewidth', 4,'color',[0 0 0 0.5]);
    text(min(SFOptions)-0.8,0.03,'No difference line','fontsize',fontSize);
    
    % Add legend.
    legendHandles = {subjectOptions{:}, 'Mean'};
    legend(legendHandles,'location','southeastoutside','fontsize',fontSize-3);
    
    % Add title.
    sgtitle(sprintf('Log sensitivity difference between Filter ( A ) and Filter ( %s )',filterOptions{ff}),'Fontsize',fontSize+5);
    
    %% Save the plot.
    if (SAVEPLOTS)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
            testFilename = fullfile(testFiledir, append(sprintf('SensitivityDifference_A_and_%s_over_SFs',filterOptions{ff}),imgFileFormat));
            saveas(gcf,testFilename);
            fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
        end
    end
end
