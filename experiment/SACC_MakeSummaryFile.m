%% SACC_MakeSummaryFile
%
% This creates a text summary file by merging CS and AUC calculation
% results for all subjects. The created file should basically contain all
% results acquired from this project.
%
% See also:
%    t_CSFGeneratorAnalyze, SACC_FitCSF

% History:
%    03/10/23   smo    Started on it.
%    04/24/23   smo    Now making CS and AUD files separately.
%    09/11/23   smo    Added sanity check for AUC file.

%% Initialize.
clear; close all;

%% Get available subjects.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFileList = dir(testFiledir);
    
    % Subject name.
    for tt = 1:length(testFileList)
        testFilenameList{tt}  = testFileList(tt).name;
    end
    idxSubjectName = find(str2double(testFilenameList)>0);
    subjectNameOptions = testFilenameList(idxSubjectName);
end

%% Read out CS text file.
nSubjects = length(subjectNameOptions);
subjectCompletedCS = {};
subjectCompletedAUC = {};
for ss = 1:nSubjects
    subjectName = subjectNameOptions{ss};
    testFilename = fullfile(testFiledir,subjectName,'CSF',sprintf('CS_Summary_%s.xlsx',subjectName));
    
    if isfile(testFilename)
        dataCS{:,ss} = readtable(testFilename);
        subjectCompletedCS{end+1} = subjectName;
    end
end

% Combine the tables.
tableCS = vertcat(dataCS{:});

% Update the numbering.
nRows = size(tableCS,1);
tableCS{:,'No'} = linspace(1,nRows,nRows)';

disp('CS data has been merged successfully!');

%% Read out AUC text file.
for ss = 1:nSubjects
    subjectName = subjectNameOptions{ss};
    testFilename = fullfile(testFiledir,subjectName,'CSF',sprintf('AUC_Summary_%s.xlsx',subjectName));
    
    if isfile(testFilename)
        dataAUC{:,ss} = readtable(testFilename);
        subjectCompletedAUC{end+1} = subjectName;
    end
end

% Combine the tables.
tableAUC = vertcat(dataAUC{:});

% Update the numbering.
nRows = size(tableAUC,1);
tableAUC{:,'No'} = linspace(1,nRows,nRows)';

disp('AUC data has been merged successfully!');

%% Sanity check for the CSF fitting data.
%
% As we found out that we used the test images using the Standard method,
% here we want to check if the CSF data was fitted within the 'good'
% contrast range of the test images.
%
% These cone contrasts are the marginal contrasts that can generate the
% desired contrast. If a test image contrast had a higher contrast than
% this, its contrast would be not perfectly modulated as desired. 
%
% These values were found using the validation code SpectralCalCheck_ver2.m.
coneContrastNormal = [-0.0427 0.0369 -0.0028];
coneContrastHigh = [-0.0579 0.0590 -0.0003];

imageContrastNormal = sqrt(sum(coneContrastNormal.^2));
imageContrastHigh = sqrt(sum(coneContrastHigh.^2));

logSensitivitySanityNormal = log10(1/imageContrastNormal);
logSensitivitySanityHigh = log10(1/imageContrastHigh);

% Read out the test image profile to figure out which test image set was
% used for 18 cpd per each subject.
testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
testFilename = GetMostRecentFileName(testFiledir,'TestImageProfile.xlsx');
imageData = readtable(testFilename);

% Find the subjects who used the high test image set in the experiment.
testImageMaxContrastHighImageSet = 0.1;
numSubjectsWithHighImageSet = imageData.Subject(find(imageData.TestImageContrastMax == testImageMaxContrastHighImageSet));

% Get the subject info.
subjectOptions = unique(tableAUC.Subject);
nSubjects = length(subjectOptions);
nSessions = size(tableAUC,1);

% Make a loop to check the sanity per each subject.
for ss = 1:nSessions
    % Get the current number of the subject.
    numSubjectTemp = table2array(tableAUC(ss,'Subject'));
    % Set the criteria differetly if high image set was used.
    if ismember(numSubjectTemp,numSubjectsWithHighImageSet)
        % High image.
        logSensitivitySanity18cpd = logSensitivitySanityHigh;
    else
        % Normal image.
        logSensitivitySanity18cpd = logSensitivitySanityNormal;
    end
    
    % Check if any sensitivity point on CSF is outside of the sanity range.
    tableVarNames = {'LogSensitivity_3cpd','LogSensitivity_6cpd','LogSensitivity_9cpd','LogSensitivity_12cpd','LogSensitivity_18cpd'};
    tableLogSensitivityTemp = table2array(tableAUC(ss,tableVarNames));
    if or(any(tableLogSensitivityTemp(1:4) < logSensitivitySanityNormal),...
            any(tableLogSensitivityTemp(end) < logSensitivitySanity18cpd))
        idxSessionBad(ss,1) = 1;
    else
        idxSessionBad(ss,1) = 0;
    end
end

% Add a new column to the table saving the findings.
tableAUC = addvars(tableAUC,idxSessionBad,'NewVariableNames','IdxSessionBad');

%% Visualize the results what we found from the above.
figure; clf; hold on;
spatialFrequencyOptions = [3 6 9 12 18];
logSpatialFrequencyOptions = log10(spatialFrequencyOptions);
spatialFrequencyOptionsStr = {'3' '6' '9' '12' '18'};

% Plot the sanity reference.
plot(log10([1 100]), [logSensitivitySanityNormal logSensitivitySanityNormal],':','linewidth',4,'color',[0 1 0 0.3]);
plot(log10([1 100]), [logSensitivitySanityHigh logSensitivitySanityHigh],':','linewidth',4,'color',[1 0 0 0.3]);

for ss = 1:nSessions
    % We will use different marker color for Normal image and High image set.
    
    % Update subject name.
    numSubjectTemp = table2array(tableAUC(ss,'Subject'));
    
    % Set marker color differently over test image set.
    markerColorNormal = 'g';
    markerColorHigh = 'r';
    if ismember(numSubjectTemp,numSubjectsWithHighImageSet)
        markerColor = markerColorHigh;
    else
        markerColor = markerColorNormal;
    end
    
    % Get CSF curve values.
    tableVarNames = {'LogSensitivity_3cpd','LogSensitivity_6cpd','LogSensitivity_9cpd','LogSensitivity_12cpd','LogSensitivity_18cpd'};
    logSensitivityTemp = table2array(tableAUC(ss,tableVarNames));
    
    % Plot it - 3, 6, 9, 12 cpd. 
    % We only used a normal image set for this spatial frequencies, so it
    % will be all same colors.
    plot(logSpatialFrequencyOptions(1:4),logSensitivityTemp(1:4),'ko','markerfacecolor',markerColorNormal);
    
    % Plot it - 18 cpd.
    % For 18 cpd, we will plot it in different colors per different image
    % set we used (normal / high).
    plot(logSpatialFrequencyOptions(end),logSensitivityTemp(end),'ko','markerfacecolor',markerColor);
end

% Set axis and stuffs.
% This is the same format as our analysis.
xlabel('Spatial Frequency (cpd)','fontsize',15);
ylabel('Log Contrast Sensitivity','fontsize',15);
xticks(logSpatialFrequencyOptions);
xticklabels(spatialFrequencyOptionsStr);
xlim([min(log10(spatialFrequencyOptions)) max(log10(spatialFrequencyOptions))]);
ylim(log10([1 600]));
yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
yticks(yaxisRange);
ytickformat('%.2f');
title(sprintf('CSF data points (N=%d)',nSessions), 'fontsize',15);
f = flip(get(gca,'children'));
legend(f([1 2 3 100]), 'Ref-Normal','Ref-High','CSF (Normal)', 'CSF (High)','location','southeast','fontsize',15);

%% Save out the result.
%
% We will create two separate excel files, one contaning the CS results and
% the other containing the AUC results.
range = 'B2';
testFilenameCS = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_CS.xlsx');
testFilenameAUC = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');

% Make CS summary file.
writetable(tableCS, testFilenameCS,'Range',range);
disp('Summary file has been saved successfully! - (CS)');

% Make AUC summary file.
writetable(tableAUC,testFilenameAUC,'Range',range);
disp('Summary file has been saved successfully! - (AUC)');
