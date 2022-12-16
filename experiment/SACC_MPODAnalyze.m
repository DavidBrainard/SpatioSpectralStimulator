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
        % Set the directory of the file saved.
        testFiledirTemp = fullfile(testFiledir,subjectNameOptions{ss});
        fileListTemp = dir(testFiledirTemp);
        
        % Read the table from the file.
        numTargetFilename = 3;
        testFilenameTemp = fullfile(testFiledirTemp,fileListTemp(numTargetFilename).name);
        data = readtable(testFilenameTemp);

        % Get the info of which eye was first measured.
        whichEyeFirstMeasure = data{3,2}{1};
        
        % Get MPOD results here.
        mpodFirstMeasure = str2double(data{6,2}{1});
        stdFirstMeasure = str2double(data{7,2}{1});
        
        % Get the info of dominant eye ('R' or 'L').
        if strcmp(subjectNameOptions(ss),'002')
            whichEyeDominant{ss} = 'L';
        else
            whichEyeDominant{ss} = data{13,2}{1};
        end
        
        % For subject 002, we are missing the data file, but we have the
        % values, so typing manually.
        if strcmp(subjectNameOptions(ss),'002')
            mpodSecondMeasure = 0.4;
            stdSecondMeasure = 0.03;
        else
            mpodSecondMeasure = str2double(data{16,2}{1});
            stdSecondMeasure = str2double(data{17,2}{1});
        end
        
        % Set the data right based on the order of measurements.
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
fontSize = 15;
xlabel('Subject name', 'fontSize', fontSize);
ylabel('MPOD', 'fontSize', fontSize);
xlim([1.5 str2double(subjectNameOptions{end})])
xticks([str2double(subjectNameOptions{1}) : 1 : str2double(subjectNameOptions{end})]);
title('MPOD results','fontsize',fontSize);
subtitle('Dominant eye is marked as triangle and big','fontsize',fontSize);

for ss = 1:nSubjectName
    numSubject = subjectNameOptions(ss);
    numSubject = str2double(numSubject);

    % For the plot, we will set different marker for dominant and
    % non-dominant eyes.
    %
    % Marker type.
    markerTypeDominantEye = '^';
    markerTypeNonDominantEye = 'o';
    
    % Marker size.
    markerSizeDominantEye = 13;
    markerSizeNonDominantEye = 5;
    
    % Marker color.
    markerColorDominantEye = 'c';
    markerColorNonDominantEye = 'k';
    
    % Set marker styles here.
    if strcmp(whichEyeDominant{ss},'L')
        markerTypeLeftEye = markerTypeDominantEye;
        markerTypeRightEye = markerTypeNonDominantEye;
        
        markerSizeLeftEye = markerSizeDominantEye;
        markerSizeRightEye = markerSizeNonDominantEye;
        
        markerColorLeftEye = markerColorDominantEye;
        markerColorRightEye = markerColorNonDominantEye;
    else
        markerTypeLeftEye = markerTypeNonDominantEye;
        markerTypeRightEye = markerTypeDominantEye;

        markerSizeLeftEye = markerSizeNonDominantEye;
        markerSizeRightEye = markerSizeDominantEye;

        markerColorLeftEye = markerColorDominantEye;
        markerColorRightEye = markerColorNonDominantEye;
    end
    
    % Plot mpod and its standard deviation.
    %
    % Left eye.
    plot(numSubject, mpodLeftEye(ss), markerTypeLeftEye,'markersize',markerSizeLeftEye,'markeredgecolor','k','markerfacecolor','b');
    errorbar(numSubject, mpodLeftEye(ss), stdLeftEye(ss), 'b','linewidth',1);

    % Right eye. 
    plot(numSubject, mpodRightEye(ss),markerTypeRightEye,'markersize',markerSizeRightEye,'markeredgecolor','k','markerfacecolor','c');
    errorbar(numSubject, mpodRightEye(ss), stdRightEye(ss), 'c','linewidth',1);

    % Connect data points between left and right eyes.
    plot([numSubject numSubject], [mpodLeftEye(ss) mpodRightEye(ss) ],'k:','LineWidth',1);
    
end
legend('Left eye (color)','','Right eye (color)', 'fontSize', fontSize)

%% Save the plot.
SAVETHEPLOT = true;
if (SAVETHEPLOT)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir,'Results_MPOD');
        testFileFormat = '.tiff';
        saveas(gcf,append(testFilename,testFileFormat));
    end
end