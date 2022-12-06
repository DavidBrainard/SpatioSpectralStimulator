% SACC_IsoLuminanceAnalyze
%
% This is for analyzing the iso-luminance determination results.
%
% See also:
%    SACC_IsoLuminanceDetermination

% History:
%    12/5/22   smo    - Started on it.

%% Initialize.
clear; close all;

%% Load the data.
if (ispref('SpatioSpectralStimulator','SACCData'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCData');
    testFileList = dir(testFiledir);

    % Find available subject name.
    for tt = 1:length(testFileList)
        testFilenameList{tt}  = testFileList(tt).name;
    end
    idxSubjectName = find(str2double(testFilenameList)>0);
    subjectNameOptions = testFilenameList(idxSubjectName);

    % Find subjects that has iso-luminance determination data.
    nSubjectName = length(subjectNameOptions);
    numSubjects = {};
    allData = {};

    for ss = 1:nSubjectName
        testFiledirTemp = fullfile(testFiledir,subjectNameOptions{ss},'FlickerPhotom');
        subjectNameTemp = subjectNameOptions{ss};

        % Load the data here if available.
        if isdir(testFiledirTemp)
            testFilenameTemp = GetMostRecentFileName(testFiledirTemp,...
                sprintf('flickerPhotom_%s',subjectNameTemp));
            data = load(testFilenameTemp);

            % Save subject name and the data.
            numSubjects{end+1} = subjectNameTemp;
            allData{end+1} = data.data;
        else
        end
    end
end

%% Plot it.
figure; clf; hold on;
xlabel('Subject name');
ylabel('Intensity matching red (8-bit)');

nSubjects = length(numSubjects);
for ss = 1:nSubjects
    numSubject = str2double(numSubjects{ss});
    dataRaw = allData{ss}.results;
    dataMean = mean(dataRaw);
    dataStd = std(dataRaw);

    % Plot raw data.
    plot(numSubject*ones(1,length(dataRaw)), dataRaw,'ko');

    % Plot mean and its standard deviation.
    plot(numSubject, dataMean,'ro','markerfacecolor','r');
    errorbar(numSubject, dataMean, dataStd, 'r');

    clear dataRaw;
end

%% Save the plot.
