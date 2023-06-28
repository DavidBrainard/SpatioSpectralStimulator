% SACC_MPODAnalyze.
%
% This is for analyzing MPOD measurement results.

% History:
%    12/05/22   smo    - Started on it.
%    02/02/23   smo    - Updated the legends on the plot.

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
        
        % Collect MPOD values of the dominant eye.
        mpodNondominant(ss) = mpodFirstMeasure;
        mpodDominant(ss) = mpodSecondMeasure;
    end
end

% Save the MPOD results.
MPOD.Subject = subjectNameOptions;
MPOD.dominantEye = mpodDominant;
MPOD.nonDominantEye = mpodNondominant;

if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = fullfile(testFiledir, 'MPOD.mat');
    save(testFilename,'MPOD');
end

%% Plot it.
figure; clf; hold on;
fontSize = 15;
xlabel('Subject name', 'fontSize', fontSize);
ylabel('MPOD', 'fontSize', fontSize);
xlim([1.5 str2double(subjectNameOptions{end})])
xticks([str2double(subjectNameOptions{1}) : 1 : str2double(subjectNameOptions{end})]);
title('MPOD results','fontsize',fontSize);
subtitle('Left eye = (triangle) / right eye = (circle)','fontsize',fontSize);

for ss = 1:nSubjectName
    % For the plot, we will set different marker for dominant and
    % non-dominant eyes.
    %
    % Marker type (left and right).
    markerTypeLeftEye = '^';
    markerTypeRightEye = 'o';
    
    % Marker size (dom or non-dom).
    markerSizeDominantEye = 7;
    markerSizeNonDominantEye = 7;
    
    % Marker color (dom or non-dom).
    markerColorDominantEye = 'r';
    markerColorNonDominantEye = 'k';
    
    % Set marker styles here.
    if strcmp(whichEyeDominant{ss},'L')
        markerSizeLeftEye = markerSizeDominantEye;
        markerSizeRightEye = markerSizeNonDominantEye;
        
        markerColorLeftEye = markerColorDominantEye;
        markerColorRightEye = markerColorNonDominantEye;
    else
        markerSizeLeftEye = markerSizeNonDominantEye;
        markerSizeRightEye = markerSizeDominantEye;
        
        markerColorLeftEye = markerColorNonDominantEye;
        markerColorRightEye = markerColorDominantEye;
    end
    
    % Plot mpod and its standard deviation.
    %
    % Left eye.
    plot(ss, mpodLeftEye(ss), markerTypeLeftEye,'markersize',markerSizeLeftEye,...
        'markeredgecolor','k','markerfacecolor',markerColorLeftEye);
    errorbar(ss, mpodLeftEye(ss), stdLeftEye(ss), markerColorLeftEye,'linewidth',1);
    
    % Right eye.
    plot(ss, mpodRightEye(ss), markerTypeRightEye,'markersize',markerSizeRightEye,...
        'markeredgecolor','k','markerfacecolor',markerColorRightEye);
    errorbar(ss, mpodRightEye(ss), stdRightEye(ss), markerColorRightEye,'linewidth',1);
    
    % Connect data points between left and right eyes.
    plot([ss ss], [mpodLeftEye(ss) mpodRightEye(ss) ],'k:','LineWidth',1);  
end
xlim([0 nSubjectName+1]);
xticks([1:1:32]);
xticklabels(subjectNameOptions);

legend('Dominant eye (red)','','Non-dominant eye (black)', 'fontSize', 13)

%% Save the plot.
SAVETHEPLOT = true;
if (SAVETHEPLOT)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir,'Results_MPOD');
        testFileFormat = '.tiff';
        saveas(gcf,append(testFilename,testFileFormat));
        disp('MPOD result plot has been saved successfully!');
    end
end