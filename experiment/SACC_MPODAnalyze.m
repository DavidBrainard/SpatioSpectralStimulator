% SACC_MPODAnalyze.
%
% This is for analyzing MPOD measurement results.

% History:
%    12/5/22   smo    - Started on it.

%% Initialize.
clear; close all;

%% Load the data.
if (ispref('SpatioSpectralStimulator','SACCMaterials'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCMaterials'),'SACC_MPOD');
    testFileList = dir(testFiledir);

    % Find available subject name.
    for tt = 1:length(testFileList)
        testFilenameList{tt}  = testFileList(tt).name;
    end
    idxSubjectName = find(str2double(testFilenameList)>0);
    subjectNameOptions = testFilenameList(idxSubjectName);

    % Read MPOD data here.
    nSubjectName = length(subjectNameOptions);

    for ss = 1:nSubjectName
        testFiledirTemp = fullfile(testFiledir,subjectNameOptions{ss});
        fileListTemp = dir(testFiledirTemp);

        numTargetFilename = 3;
        testFilenameTemp = fullfile(testFiledirTemp,fileListTemp(numTargetFilename).name);
        data = readtable(testFilenameTemp);

        % Get the data of both eyes.
        whichEyeFirstMeasure = data{3,2}{1};

        mpodFirstMeasure = str2double(data{6,2}{1});
        stdFirstMeasure = str2double(data{7,2}{1});

        if strcmp(subjectNameOptions(ss),'002')
            mpodSecondMeasure = 0.4;
            stdSecondMeasure = 0.03;
        else
            mpodSecondMeasure = str2double(data{16,2}{1});
            stdSecondMeasure = str2double(data{17,2}{1});
        end

        if (strcmp(whichEyeFirstMeasure,'L'))
            mpodLeftEye(ss) = mpodFirstMeasure;
            stdLeftEye(ss) = stdFirstMeasure;
            mpodRightEye(ss) = mpodSecondMeasure;
            stdRightEye(ss) = stdSecondMeasure;

        elseif (strcmp(whichEyeFirstMeasure,'R'))
            mpodLeftEye(ss) = mpodSecondMeasure;
            stdLeftEye(ss) = stdSecondMeasure;
            mpodRightEye(ss) = mpodFirstMeasure;
            stdRightEye(ss) = stdFirstMeasure;
        end
    end
end

%% Plot it.
figure; clf; hold on;
xlabel('Subject name');
ylabel('MPOD');

for ss = 1:nSubjectName
    numSubject = subjectNameOptions(ss);
    numSubject = str2double(numSubject);

    % Plot mpod and its standard deviation.
    %
    % Left eye.
    plot(numSubject, mpodLeftEye(ss),'ro','markerfacecolor','r');
    errorbar(numSubject, mpodLeftEye(ss), stdLeftEye(ss), 'r');

    % Right eye.
    plot(numSubject, mpodRightEye(ss),'go','markerfacecolor','g');
    errorbar(numSubject, mpodRightEye(ss), stdRightEye(ss), 'g');

    % Connect data points between left and right eyes.
    plot([numSubject numSubject], [mpodLeftEye(ss) mpodRightEye(ss) ],'k--','LineWidth',1);
    
end
legend('Left eye','','Right eye')

%% Save the plot.
