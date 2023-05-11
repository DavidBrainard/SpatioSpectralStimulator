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

%% Set the options. 
SAVETHEPLOT = true;

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
    allCalData = {};
    
    for ss = 1:nSubjectName
        testFiledirTemp = fullfile(testFiledir,subjectNameOptions{ss},'FlickerPhotom');
        subjectNameTemp = subjectNameOptions{ss};
        
        % Load the data here if available.
        if isdir(testFiledirTemp)
            testFilenameTemp = GetMostRecentFileName(testFiledirTemp,...
                sprintf('flickerPhotom_%s',subjectNameTemp));
            data = load(testFilenameTemp);
            
            % Load the calibration file used for the experiment.
            [dir filename format] = fileparts(testFilenameTemp);
            numExtract = regexp(filename,'\d+','match');
            date = sprintf('%s-%s-%s',numExtract{2},numExtract{3},numExtract{4});
            
            calFiledir = fullfile(testFiledir,'CheckFlickerPhotom');
            calDataFilename = GetMostRecentFileName(calFiledir,...
                sprintf('checkFlickerPhotom_%s',date));
            calData = load(calDataFilename);
        
            % Save subject name and the data.
            numSubjects{end+1} = subjectNameTemp;
            allData{end+1} = data.data;
            allCalData{end+1} = calData;
        end
    end
end

%% Plot the raw data.
figure; clf; hold on;
fontSize = 15;
xlabel('Subject name','fontsize',fontSize);
ylabel('Output (no unit)','fontsize',fontSize);
xlim([1.5 str2double(numSubjects{end})])
ylim([0 1.4]);
xticks([str2double(numSubjects{1}) : 1 : str2double(numSubjects{end})]);
title('Iso-luminance determination results','fontsize',fontSize);
subtitle('This is the raw results');

nSubjects = length(numSubjects);
for ss = 1:nSubjects
    numSubject = str2double(numSubjects{ss});
    
    if (numSubject == 2)
        redInputSettings = [92 104 116 104 100 112];
    else
        redInputSettings = allData{ss}.results;
    end
    
    % Linear fitting to get output settings of the matching results (no unit).
    x = allCalData{ss}.data.redIntensity;
    y = allCalData{ss}.dataBC.xyzRedBC(2,:);
    p = polyfit(x, y, 1);
    greenOutputSetting = allCalData{ss}.dataBC.xyzGreenBC(2);
    redOutputSettings = polyval(p, redInputSettings);
    
    % Calculate mean and std of the output settings.
    meanRedOutputSettings = mean(redOutputSettings);
    stdRedOutputSettings = std(redOutputSettings);
    
    % Plot raw data (Red).
    plot(ss*ones(1,length(redOutputSettings)), redOutputSettings,'ko');
    
    % Plot the reference (Green).
    plot(ss, greenOutputSetting,'go','markersize',9,'markeredgecolor','k','markerfacecolor','g');
    
    % Plot mean and its standard deviation (Red).
    plot(ss, meanRedOutputSettings,'ro','markersize',9,'markeredgecolor','k','markerfacecolor','r');
    errorbar(ss, meanRedOutputSettings, stdRedOutputSettings, 'r', 'linewidth', 1);
end
xlim([0,nSubjects+1]);
xticks([1:1:nSubjects]);
xticklabels(numSubjects);
legend('Raw data (red)','Reference (green)','Mean (red)','fontsize',fontSize,'location','southeast');

% Save the plot.
if (SAVETHEPLOT)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir,'Results_FlickerPhotom_raw');
        testFileFormat = '.tiff';
        saveas(gcf,append(testFilename,testFileFormat));
        disp('Iso-luminance determination result plot has been successfully saved! - (raw)');
    end
end

%% Plot the normalized data.
figure; clf; hold on;
fontSize = 15;
xlabel('Subject name','fontsize',fontSize);
ylabel('Output (no unit)','fontsize',fontSize);
xlim([1.5 str2double(numSubjects{end})])
ylim([0 1.4]);
xticks([str2double(numSubjects{1}) : 1 : str2double(numSubjects{end})]);
title('Iso-luminance determination results','fontsize',fontSize);
subtitle('This is the normalized results for reference(green) to be 1');

nSubjects = length(numSubjects);
for ss = 1:nSubjects
    numSubject = str2double(numSubjects{ss});
    
    if (numSubject == 2)
        redInputSettings = [92 104 116 104 100 112];
    else
        redInputSettings = allData{ss}.results;
    end
    
    % Linear fitting to get output settings of the matching results (no unit).
    x = allCalData{ss}.data.redIntensity;
    y = allCalData{ss}.dataBC.xyzRedBC(2,:);
    p = polyfit(x, y, 1);
    greenOutputSetting = allCalData{ss}.dataBC.xyzGreenBC(2);
    redOutputSettings = polyval(p, redInputSettings);
    
    % Normalize it if you want.
    greenOutputSetting = greenOutputSetting/greenOutputSetting;
    redOutputSettings = redOutputSettings./greenOutputSetting;
    
    % Calculate mean and std of the output settings.
    meanRedOutputSettings = mean(redOutputSettings);
    stdRedOutputSettings = std(redOutputSettings);
    
    % Plot raw data (Red).
    plot(ss*ones(1,length(redOutputSettings)), redOutputSettings,'ko');
    
    % Plot the reference (Green).
    plot(ss, greenOutputSetting,'go','markersize',9,'markeredgecolor','k','markerfacecolor','g');
    
    % Plot mean and its standard deviation (Red).
    plot(ss, meanRedOutputSettings,'ro','markersize',9,'markeredgecolor','k','markerfacecolor','r');
    errorbar(ss, meanRedOutputSettings, stdRedOutputSettings, 'r', 'linewidth', 1);
end
xlim([0,nSubjects+1]);
xticks([1:1:nSubjects]);
xticklabels(numSubjects);
legend('Raw data (red)','Reference (green)','Mean (red)','fontsize',fontSize,'location','southeast');

% Save the plot.
if (SAVETHEPLOT)
    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = fullfile(testFiledir,'Results_FlickerPhotom_normalized');
        testFileFormat = '.tiff';
        saveas(gcf,append(testFilename,testFileFormat));
        disp('Iso-luminance determination result plot has been successfully saved! - (normalized)');
    end
end