% SACC_FitCSFReadOnly
%
% This is to read the existing fitting results.
%
% See also:
%    SACC_FitCSFReadOnly.

% History:
%    3/28/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
fittingMode = 'crossValBootAcrossFmincon';

%% Check the subjects with all data.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFileList = dir(testFiledir);
    
    % Find available data of subject and spatial frequency.
    for tt = 1:length(testFileList)
        testFilenameList{tt}  = testFileList(tt).name;
    end
    idxSubjectName = find(str2double(testFilenameList)>0);
    subjectNameOptions = testFilenameList(idxSubjectName);
    
    % Get the subjects having all figures available.
    subjectAvailable = {};
    for ss = 1:length(subjectNameOptions)
        % Get the subject number.
        whichSub = subjectNameOptions{ss};
        
        % Set the directory of the figures saved per each subject.
        testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
        testFiledir = fullfile(testFiledir,whichSub,'CSF');
        targetFilename = append(sprintf('CSF_%s',whichSub));
        fileList = dir(append(fullfile(testFiledir,targetFilename),'*'));
        
        % Check the number of available figures per subject.
        nFiles = 0;
        for ff = 1:length(fileList)
            fileNameTemp = string(fileList(ff).name);
            nFiles = nFiles + contains(fileNameTemp,fittingMode);
        end
        
        % Save if the subject has all available figures.
        if (nFiles == 5)
            subjectAvailable{end+1} = whichSub;
        end
    end
end

%% Prompt which subject result to read.
while 1
    for 
    textSubAvailable = 
    end 
    inputMessage = sprintf('Which subject to test: \n %s',);
    whichSub = input(inputMessage, 's');
    whichFilterOptions = {'A', 'B', 'C', 'D', 'E'};
    
    if ismember(whichSub, whichFilterOptions)
        break
    end
    
    disp('Filter should be chose within [A, B, C, D, E]!');
end

%% Display the result.
