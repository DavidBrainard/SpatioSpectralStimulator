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
% for ss = 1:nSFs
%     logSensitivityAllFiltersTemp = squeeze(logSensitivity(:,ss,:));
%     logSensitivityAllFiltersTemp = reshape(logSensitivityAllFiltersTemp, [size(logSensitivityAllFiltersTemp,1)*size(logSensitivityAllFiltersTemp,2) 1]);
%     stdErrorAllFiltersAndSubjects(ss) = std(logSensitivityAllFiltersTemp)/sqrt(length(logSensitivityAllFiltersTemp));
% end

%% Get subject gender info, MPOD, Iso-luminance results.
%
% Subject gender.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'Experiment_Progress.xlsx');
    tableGender = readtable(testFilename,'UseExcel',true,'sheet','SACC_Exp_Progress','range','B2:C48');
end

% Find the index for male and female out of all subjects.
for ss = 1:nSubjects
    % Read out subject gender from the table.
    genderSubject = tableGender.SexAAB(tableGender.SubjectNo_(str2num(subjectOptions{ss})));
    genderSubject = genderSubject{:};
    
    % Allocate index number according to gender.
    male = 1;
    female = 0;
    switch genderSubject
        case 'M'
            isMale = male;
        case 'F'
            isMale = female;
    end
    
    % Save out index for male.
    genderSubjects(ss) = isMale;
end
% Convert the index into an array with the subject order.
idxSubjectMale = find(genderSubjects == male);
idxSubjectFemale = find(genderSubjects == female);

% MPOD.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'MPOD.mat');
    dataMPOD = load(testFilename);
end

% Read out dominant eye MPOD data except subject 021 as it is weirdly high.
MPOD = dataMPOD.MPOD.dominantEye;

% Find the index of subject 021.
for ss = 1:nSubjects
    numSubjectTemp = dataMPOD.MPOD.Subject{ss};
    if strcmp(numSubjectTemp,'021')
        idxMPODSwitch = ss;
    end
end

% Switch MPOD value of subject 021 to non-dominant result.
MPOD(idxMPODSwitch) = dataMPOD.MPOD.nonDominantEye(idxMPODSwitch);

% Iso-luminance determination results.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'IsoLuminanceDetermination.mat');
    dataFlicker = load(testFilename);
end

% Read out data we need.
meanRedNormalizedFlicker = dataFlicker.Flicker.MeanRedNormalized;
meanGreenFlicker = mean(cell2mat(dataFlicker.Flicker.RawDataGreen));

% Calculate flicker score. This will be used for the analysis.
flickerScore = meanRedNormalizedFlicker./meanGreenFlicker;

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

%% 4-1) Plot Grand cCSF (line graph - separate).
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

%% 4-2) Grand grand mean cCSF.
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

%% 4-3) Plot Grand CCSF (line graph - together in one figure).
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
    
    % Set the figure and font size on axes.
    figPosition = [0 0 1200 1000];
    fontSize = 15;
    
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
    
    % Standard Error.
    errorbar(SFOptions, mean(logSensitivityFilterTest-logSensitivityFilterA), std(logSensitivityFilterTest-logSensitivityFilterA)/sqrt(length(logSensitivityFilterTest)),...
        'b','linewidth',2);
    ee = get(gca,'Children');
    set(ee(1),'LineStyle','none');
    
    % No difference line.
    plot([min(SFOptions)-1 max(SFOptions+1)], [0 0], 'k-','linewidth', 4,'color',[0 0 0 0.5]);
    text(min(SFOptions)-0.8,0.05,'No difference line','fontsize',fontSize);
    
    % Add legend.
    legendHandles = {subjectOptions{:}, 'Mean'};
    legend(legendHandles,'location','southeastoutside','fontsize',8);
    
    % Add title.
    sgtitle(sprintf('Log sensitivity difference between Filter ( A ) and Filter ( %s )',filterOptions{ff}),'Fontsize',fontSize+2);
    
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

%% 6-1) Gender effect on cCSF.
% Same mean CSF graph, but making an average separately of male and female.
%
% Separate male and female sensitiivty data.
logSensitivityMale = logSensitivity(idxSubjectMale,:,:);
logSensitivityFemale = logSensitivity(idxSubjectFemale,:,:);

% Make an average over the filters and subjects.
% Male.
for ss = 1:nSFs
    meanLogSensitivityMaleTemp = squeeze(logSensitivityMale(:,ss,:));
    meanLogSensitivityMaleTemp = reshape(meanLogSensitivityMaleTemp,size(meanLogSensitivityMaleTemp,2)*size(meanLogSensitivityMaleTemp,1),1);
    meanLogSensitivityMalePerSFs(ss) = mean(meanLogSensitivityMaleTemp);
    stdErrorMaleOverFiltersPerSFs(ss) = std(meanLogSensitivityMaleTemp)/sqrt(length(meanLogSensitivityMaleTemp));
end

% Female.
for ss = 1:nSFs
    meanLogSensitivityFemaleTemp = squeeze(logSensitivityFemale(:,ss,:));
    meanLogSensitivityFemaleTemp = reshape(meanLogSensitivityFemaleTemp,size(meanLogSensitivityFemaleTemp,2)*size(meanLogSensitivityFemaleTemp,1),1);
    meanLogSensitivityFemalePerSFs(ss) = mean(meanLogSensitivityFemaleTemp);
    stdErrorFemaleOverFiltersPerSFs(ss) = std(meanLogSensitivityFemaleTemp)/sqrt(length(meanLogSensitivityFemaleTemp));
end

% Plot it.
figure; hold on;
plot(SFOptions, meanLogSensitivityMalePerSFs,'b.-','color',[0 0 1 0.5],'linewidth',2);
plot(SFOptions, meanLogSensitivityFemalePerSFs,'r.-','color',[1 0 0 0.5],'linewidth',2);
errorbar(SFOptions, meanLogSensitivityMalePerSFs,stdErrorMaleOverFiltersPerSFs,'color',[0 0 1 0.5]);
errorbar(SFOptions, meanLogSensitivityFemalePerSFs,stdErrorFemaleOverFiltersPerSFs,'color',[1 0 0 0.5]);
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity (N=%d)', length(meanLogSensitivityMaleTemp)),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
nMales = length(idxSubjectMale);
nFemales = length(idxSubjectFemale);
legend(sprintf('Male (N=%d)',nMales),sprintf('Female (N=%d)',nFemales),'fontsize',15);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('GenderEffect_CCSF',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 6-2) Gender effect on AUC.
% Separate the AUC data between male and female.
AUCMale = AUC(idxSubjectMale,:);
AUCFemale = AUC(idxSubjectFemale,:);
meanAUCMale = mean(AUCMale,1);
meanAUCFemale = mean(AUCFemale,1);

figure; clf; hold on;
figPosition = [0 0 900 700];
set(gcf,'position',figPosition);

xAxisTicks = [1:nFilters];
xticks(xAxisTicks);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('AUC','Fontsize',15);

% Plot male individual results.
for ss = 1:nMales
    plot(xAxisTicks,AUCMale(ss,:),'bo-','linewidth',0.7);
end

% Plot female individual results.
for ss = 1:nFemales
    plot(xAxisTicks,AUCFemale(ss,:),'ro-','linewidth',0.7);
end

% Plot mean results.
plot(xAxisTicks,meanAUCMale,'k+-','color',[0 0 1 0.5],'linewidth',7);
plot(xAxisTicks,meanAUCFemale,'k+-','color',[1 0 0 0.5],'linewidth',7);

% Add texts for the mean values.
text(xAxisTicks, meanAUCMale+0.3, num2str(round(meanAUCMale,2)'),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'fontsize', 11, 'fontweight', 'bold', 'color', 'k');
text(xAxisTicks, meanAUCFemale+0.3, num2str(round(meanAUCFemale,2)'),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'fontsize', 11, 'fontweight', 'bold', 'color', 'k');

% Some plot formats.
title('Mean AUC results: Male vs. Female','Fontsize',15);
subtitle(sprintf('Mean result is the average of (%d) subjects per gender',nMales),'FontSize',13);
f = flip(get(gca, 'children'));

legendHandles = {'Male (Individual)', 'Female (Individual)','Mean (Male)', 'Mean (Female)'};
legend(f([1 17 33 34]),legendHandles,'location','northeastoutside','fontsize',14);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('GenderEffect_AUC',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 7-1) MPOD on CCSF - Group results.
%
% We will make three separate subject groups based on the MPOD values, low,
% medium, and high based on the percentile of 33% and 66%.
MPOD_33 = prctile(MPOD,33);
MPOD_66 = prctile(MPOD,66);
MPOD_low = MPOD(MPOD <= MPOD_33);
MPOD_med = MPOD(and(MPOD > MPOD_33,MPOD <= MPOD_66));
MPOD_high = MPOD(MPOD > MPOD_66);

meanMPOD_low = mean(MPOD_low);
meanMPOD_med = mean(MPOD_med);
meanMPOD_high = mean(MPOD_high);

% Check if the numbers are right.
MPOD_recollect = sort([MPOD_low MPOD_med MPOD_high],'ascend');
MPOD_check = sort(MPOD,'ascend');
if any(~MPOD_recollect == MPOD_check)
    error('MPOD number is mismatched');
end

% Get the index.
idxMPOD_low = find(MPOD <= MPOD_33);
idxMPOD_med = find(and(MPOD > MPOD_33, MPOD <= MPOD_66));
idxMPOD_high = find(MPOD > MPOD_66);

% Get number of elements per group.
nLows = length(idxMPOD_low);
nMeds = length(idxMPOD_med);
nHighs = length(idxMPOD_high);

% Separate the log sensitivity over the MPOD groups.
logSensitivity_MPOD_low = logSensitivity(idxMPOD_low,:,:);
logSensitivity_MPOD_med = logSensitivity(idxMPOD_med,:,:);
logSensitivity_MPOD_high = logSensitivity(idxMPOD_high,:,:);

% Make an average over the group and subjects.
% MPOD low.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_MPOD_low(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_MPOD_low(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_MPOD_low(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end

% MPOD medium.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_MPOD_med(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_MPOD_med(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_MPOD_med(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end

% MPOD high.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_MPOD_high(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_MPOD_high(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_MPOD_high(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end

% Plot it.
figure; hold on;
plot(SFOptions, meanLogSensitiivty_MPOD_low,'k.-','color',[1 0 0 0.5],'linewidth',2);
plot(SFOptions, meanLogSensitiivty_MPOD_med,'k.-','color',[0.1 0.7 0.1 0.5],'linewidth',2);
plot(SFOptions, meanLogSensitiivty_MPOD_high,'k.-','color',[0 0 1 0.5],'linewidth',2);
errorbar(SFOptions, meanLogSensitiivty_MPOD_low,stdErrorLogSensitivity_MPOD_low,'color',[1 0 0 0.5]);
errorbar(SFOptions, meanLogSensitiivty_MPOD_med,stdErrorLogSensitivity_MPOD_med,'color',[0.1 0.7 0.1 0.5]);
errorbar(SFOptions, meanLogSensitiivty_MPOD_high,stdErrorLogSensitivity_MPOD_high,'color',[0 0 1 0.5]);
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity'),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
legend(sprintf('MPOD low %.2f (N=%d)',meanMPOD_low,nLows),sprintf('MPOD med %.2f (N=%d)',meanMPOD_med,nMeds),...
    sprintf('MPOD high %.2f (N=%d)',meanMPOD_high,nHighs),'fontsize',15);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('MPOD_CCSF_Group',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 7-2) MPOD on CCSF - Individual results.
meanLogSensitiivtyAcrossFilters = mean(logSensitivity,3);
figure; hold on;
MPODSorted = unique(sort(MPOD,'ascend'));
lineColors = linspace(0,1,length(MPODSorted));

for ss = 1:nSubjects
    idxMPODSorted = find(MPODSorted == MPOD(ss));
    lineColorTemp = lineColors(idxMPODSorted);
    plot(SFOptions, meanLogSensitiivtyAcrossFilters(ss,:),'k.-','color',[lineColorTemp lineColorTemp 0 0.5],'linewidth',2);
end
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity'),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
title('MPOD vs. Individual CCSF','fontsize',15);
subtitle('More yellowish line is, higher MPOD value', 'fontsize',15);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('MPOD_CCSF_Individual',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 7-3) MPOD effect on AUC.
%
% MPOD vs. AUC over the filters.
figure;
figPosition = [0 0 1500 900];
set(gcf,'position',figPosition);
for ff = 1:nFilters
    subplot(2,3,ff); hold on;
    plot(MPOD, AUC(:,ff), 'ko', 'markerfacecolor', [0.2*ff 0.2*ff 0],'markersize',8);
    xlabel('MPOD','fontsize',15);
    ylabel('AUC','fontsize',15);
    ylim([20 36]);
    title(append('Filter ',filterOptions{ff}), 'fontsize',15);
end

% Calculate the mean AUC across filters.
meanAUCacrossFilters = mean(AUC,2);

% Standard error per each point.
for ss = 1:nSubjects
    stdErrorAUCacrossFilters(ss) = std(AUC(ss,:))./sqrt(length(AUC(ss,:)));
end

% MPOD vs. mean AUC across the filters.
subplot(2,3,6); hold on;
plot(MPOD,meanAUCacrossFilters,'ro','markeredgecolor','k','markerfacecolor',[0 0 1],'markersize',8);
errorbar(MPOD,meanAUCacrossFilters,stdErrorAUCacrossFilters,'b');
a = get(gca,'children');
a(1).LineStyle = 'none';
xlabel('MPOD','fontsize',15);
ylabel('Mean AUC','fontsize',15);
ylim([20 36]);
legend('Mean','fontsize',15);
title('Mean AUC across filters', 'fontsize',15);
subtitle('Along with standard error bar','fontsize',15);

% Big title of the figure.
sgtitle('MPOD vs. AUC', 'fontsize',18);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('MPOD_AUC',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 8-1) Iso-luminance determination results on CCSF - Group results.
%
% Get index for separating groups.
idxFlickerLow = find(flickerScore <= prctile(flickerScore,33));
idxFlickerMed = find(and(flickerScore > prctile(flickerScore,33), flickerScore <= prctile(flickerScore,66)));
idxFlickerHigh = find(flickerScore > prctile(flickerScore,66));

% Separate the sensitivity by groups.
logSensitivity_flicker_low = logSensitivity(idxFlickerLow,:,:);
logSensitivity_flicker_med = logSensitivity(idxFlickerMed,:,:);
logSensitivity_flicker_high = logSensitivity(idxFlickerHigh,:,:);

% Make an average over the group and subjects.
% Flicker low.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_flicker_low(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_flicker_low(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_flicker_low(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end

% Flicker medium.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_flicker_med(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_flicker_med(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_flicker_med(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end

% Flicker high.
for ss = 1:nSFs
    logSensitivityTemp = squeeze(logSensitivity_flicker_high(:,ss,:));
    logSensitivityTemp = reshape(logSensitivityTemp,size(logSensitivityTemp,2)*size(logSensitivityTemp,1),1);
    meanLogSensitiivty_flicker_high(ss) = mean(logSensitivityTemp);
    stdErrorLogSensitivity_flicker_high(ss) = std(logSensitivityTemp)/sqrt(length(logSensitivityTemp));
end


% Plot it.
figure; hold on;
plot(SFOptions, meanLogSensitiivty_flicker_low,'k.-','color',[1 0 0 0.5],'linewidth',2);
plot(SFOptions, meanLogSensitiivty_flicker_med,'k.-','color',[0.1 0.7 0.1 0.5],'linewidth',2);
plot(SFOptions, meanLogSensitiivty_flicker_high,'k.-','color',[0 0 1 0.5],'linewidth',2);
errorbar(SFOptions, meanLogSensitiivty_flicker_low,stdErrorLogSensitivity_flicker_low,'color',[1 0 0 0.5]);
errorbar(SFOptions, meanLogSensitiivty_flicker_med,stdErrorLogSensitivity_flicker_med,'color',[0.1 0.7 0.1 0.5]);
errorbar(SFOptions, meanLogSensitiivty_flicker_high,stdErrorLogSensitivity_flicker_high,'color',[0 0 1 0.5]);
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity'),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
legend(sprintf('Flicker low %.2f (N=%d)',mean(flickerScore(idxFlickerLow)),length(idxFlickerLow)),...
    sprintf('Flicker med %.2f (N=%d)',mean(flickerScore(idxFlickerMed)),length(idxFlickerMed)),...
    sprintf('Flicker high %.2f (N=%d)',mean(flickerScore(idxFlickerHigh)),length(idxFlickerHigh)),...
    'fontsize',15);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('IsoLuminace_CCSF_Group',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 8-2) Iso-luminance determination results on CCSF - Individual.
meanLogSensitiivtyAcrossFilters = mean(logSensitivity,3);
figure; hold on;
flickerScoreSorted = unique(sort(flickerScore,'ascend'));
lineColors = linspace(0,1,length(flickerScoreSorted));

for ss = 1:nSubjects
    idxFlickerScoreSorted = find(flickerScoreSorted == flickerScore(ss));
    lineColorTemp = lineColors(idxFlickerScoreSorted);
    plot(SFOptions, meanLogSensitiivtyAcrossFilters(ss,:),'k.-','color',[lineColorTemp 0 0 0.5],'linewidth',2);
end
xticks(SFOptions);
xticklabels(SFOptions);
xlabel('Spatial Frequency','Fontsize',15);
ylabel(sprintf('Mean Log Contrast Sensitivity'),'Fontsize',15);
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
ylim(log10([1 600]));
yticks(yaxisRange);
ytickformat('%.2f');
title('Flicker score vs. Individual CCSF','fontsize',15);
subtitle('More reddish line is, higher flicker score', 'fontsize',15);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('IsoLuminance_CCSF_Individual',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end

%% 8-3) Iso-luminance determination results on AUC - Individual.
%
% Iso-luminance vs. AUC over the filters.
figure;
figPosition = [0 0 1500 900];
set(gcf,'position',figPosition);
for ff = 1:nFilters
    subplot(2,3,ff); hold on;
    plot(flickerScore, AUC(:,ff), 'ko', 'markerfacecolor', [0.2*ff 0.2*ff 0],'markersize',8);
    plot([1 1], [0 40], 'k-','color', [0.5 0.5 0.5 0.2], 'linewidth',5);
    xlabel('Flicker score','fontsize',15);
    ylabel('AUC','fontsize',15);
    ylim([20 36]);
    title(append('Filter ',filterOptions{ff}), 'fontsize',15);
end

% Calculate the mean AUC across filters.
meanAUCacrossFilters = mean(AUC,2);

% Standard error per each point.
for ss = 1:nSubjects
    stdErrorAUCacrossFilters(ss) = std(AUC(ss,:))./sqrt(length(AUC(ss,:)));
end

% Flicker score vs. mean AUC across the filters.
subplot(2,3,6); hold on;
plot(flickerScore,meanAUCacrossFilters,'ro','markeredgecolor','k','markerfacecolor',[1 0 0],'markersize',8);
plot([1 1], [0 40], 'k-','color', [0.5 0.5 0.5 0.2], 'linewidth',5);
errorbar(flickerScore,meanAUCacrossFilters,stdErrorAUCacrossFilters,'r');
a = get(gca,'children');
a(1).LineStyle = 'none';
xlabel('Flicker Score','fontsize',15);
ylabel('Mean AUC','fontsize',15);
ylim([20 36]);
legend('Mean','fontsize',15);
title('Mean AUC across filters', 'fontsize',15);
subtitle('Along with standard error bar','fontsize',15);

% Big title of the figure.
sgtitle('Iso-luminance determination vs. AUC', 'fontsize',18);

% Save the plot.
if (SAVEPLOTS)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir, append('IsoLuminance_AUC',imgFileFormat));
        saveas(gcf,testFilename);
        fprintf('Plot has been saved successfully! - (%s) \n', testFilename);
    end
end
