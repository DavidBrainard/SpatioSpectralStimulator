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
% Bar graph.
figure; hold on;
xAxisBar = [1:length(meanAUC)];
bar(xAxisBar,meanAUC);
xticks(xAxisBar);
xticklabels(filterOptions);
xlabel('Filter','Fontsize',15);
ylabel('Mean AUC','Fontsize',15);
errorbar(xAxisBar,meanAUC,stdAUC,'k','linewidth',1);
legend(sprintf('Mean AUC (N=%d)',nAUCPerSub),'FontSize',13);
title('Mean AUC over the Filters','Fontsize',15);
