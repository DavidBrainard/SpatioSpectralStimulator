% SACC_AUCAnalyze
%
% This is to make a summary plot of AUC of all subjects according to
% different filters.

% History:
%    5/2/23    smo    - Wrote it.

%% Initialize.
clear; close all;

%% Read the AUC summary table.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
    T = readtable(testFilename);
end

%% Sort the AUC data over the filters.
%
% Sort the table over the filter type.
sortingCriteria = 'Filter';
T = sortrows(T,sortingCriteria);

% Filter options.
filterOptions = {'A', 'B', 'C', 'D', 'E'};
nFilters = length(filterOptions);

% Subjects.
nSubjects = length(unique(T.Subject));

% Count the number of subjects per each filter. This part will be
% elaborated later on. For now, it works fine.
nAUCPerSub = 0;
for nn = 1:size(T,1)
    if (T.Filter{nn} == filterOptions{1})
        nAUCPerSub = nAUCPerSub + 1;
    end
end

% Extract AUC over the filters.
for ff = 1:nFilters
    AUCPerFilter(:,ff) = T.AUC(1+nAUCPerSub*(ff-1):nAUCPerSub*ff);
end

% Calculate mean and standard deviations.
meanAUC = mean(AUCPerFilter,1);
stdAUC = std(AUCPerFilter,1);

%% Plot the results.
%
% 1) Mean AUC results - bar graph.
figure; hold on;
xAxisBar = [1:length(meanAUC)];
bar(xAxisBar,meanAUC);
xticks(xAxisBar);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('Mean AUC','Fontsize',15);
errorbar(xAxisBar,meanAUC,stdAUC,'k','linewidth',1);
legend(sprintf('Mean AUC (N=%d)',nAUCPerSub),'FontSize',13);
title('Mean AUC results over the Filters','Fontsize',15);

% 2) AUC over the subjects - bar graph.
indvFig = figure; clf; hold on;
figPosition = [100 100 1500 1000];
set(indvFig, 'position', figPosition);
subjectNames = T.Subject(1:nAUCPerSub);
for ff = 1:nFilters
    subplot(2,3,ff);
    xAxisBar = [1:size(AUCPerFilter,1)];
    bar(xAxisBar,AUCPerFilter(:,ff));
    xticks(xAxisBar);
    xticklabels(subjectNames);
    xlabel('Subject','Fontsize',15);
    ylabel('AUC','Fontsize',15);
    title(sprintf('Filter = %s',filterOptions{ff}),'fontsize',15);
end
sgtitle('AUC results of individual subjects per each Filter','Fontsize',20);

% 3) AUC results over the filters - line graph.
figure; clf; hold on;
xAxisBar = [1:nFilters];
xticks(xAxisBar);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('AUC','Fontsize',15);
for ss = 1:nSubjects
    plot(xAxisBar,AUCPerFilter(ss,:),'o-','linewidth',0.7);
end
plot(xAxisBar,meanAUC,'k+-','color',[0 0 0 0.4],'linewidth',7);
legendLinePlot = append('Sub ',subjectNames);
legendLinePlot{end+1} = 'Mean';
legend(legendLinePlot,'location','northeastoutside');
title('AUC results over the filters per each subject','fontsize',15);
