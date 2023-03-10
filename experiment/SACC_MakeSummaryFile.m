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

%% Initialize.
clear; close all;

% %% Find available subjects.
% %
% % Here we find out the subjects who completed study by reading out the
% % Matlab summary file.
% if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
%     testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
%     testFilename = fullfile(testFiledir,'CSFAnalysisOutput');
%     theData = load(testFilename);
%
%     % Close the plots if any pops up.
%     close all;
% else
%     error('Cannot find the data file!');
% end
%
% % Get some useful info.
% subjectName = theData.subjectNameOptions;
% sineFreqCyclesPerDeg = theData.spatialFrequencyOptions;
%
% % Find subjects who completed the study.
% nSubjects = size(sineFreqCyclesPerDeg,2);
% maxNSpatialFrequencies = theData.maxNSpatialFrequencies;
% for ss = 1:nSubjects
%     for dd = 1:maxNSpatialFrequencies
%         nSFCompleted(dd,ss) = ~isempty(sineFreqCyclesPerDeg{dd,ss});
%     end
% end
% nSFCompleted = sum(nSFCompleted,1);
% index = find(nSFCompleted==maxNSpatialFrequencies);
%
% % Take only subjects with all data.
% subjectName = subjectName(index);
% nSubjects = size(subjectName,2);

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
% We will create one new excel file that contains all data loaded here with
% two separate spreadsheets, one with CS results and the other containing
% the AUC results.
% Write a table to the excel file.
sheetCS = 1;
sheetAUC = 2;
range = 'B2';
testFilename = fullfile(testFiledir, 'SACC_Experiment_Results_Summary.xlsx');

writetable(tableCS, testFilename,'Sheet','CS', 'Range',range);
writetable(tableAUC,testFilename,'Sheet','AUC','Range',range);

disp('Summary (CS+AUC) file has been saved successfully!');
