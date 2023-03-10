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
for ss = 1:nSubjects
    subjectName = subjectNameOptions{ss};
    testFilename = fullfile(testFiledir,subjectName,'CSF',sprintf('CS_Summary_%s.xlsx',subjectName));
      
    if isfile(testFilename)
       dataCS = load(testFilename);
    end
end

%% Read out AUC text file.


%% Save out the result.

