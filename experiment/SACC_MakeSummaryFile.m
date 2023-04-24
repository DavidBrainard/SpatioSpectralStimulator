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
