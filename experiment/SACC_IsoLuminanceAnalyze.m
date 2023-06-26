% SACC_IsoLuminanceAnalyze
%
% This is for analyzing the iso-luminance determination results.
%
% See also:
%    SACC_IsoLuminanceDetermination, LightConversions, SPDToLuminance.

% History:
%    12/05/22   smo    - Started on it.
%    06/07/23   smp    - Now we plot the results using the luminance value.

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
ylabel('Luminance (cd/m2)','fontsize',fontSize);
xlim([1.5 str2double(numSubjects{end})])
% ylim([0 1.4]);
xticks([str2double(numSubjects{1}) : 1 : str2double(numSubjects{end})]);
title('Iso-luminance determination results','fontsize',fontSize);
subtitle('Raw results');

nSubjects = length(numSubjects);
for ss = 1:nSubjects
    numSubject = str2double(numSubjects{ss});
    
    if (numSubject == 2)
        redInputSettings = [92 104 116 104 100 112];
    else
        redInputSettings = allData{ss}.results;
    end
    
% We used to do use arb power with no unit, but now we have updated to use
% the luminance value as below.
%     % Linear fitting to get output settings of the matching results (no unit).
%     x = allCalData{ss}.data.redIntensity;
%     y = allCalData{ss}.dataBC.xyzRedBC(2,:);
%     p = polyfit(x, y, 1);
%     redOutputSettings = polyval(p, redInputSettings);
%     greenOutputSetting = allCalData{ss}.dataBC.xyzGreenBC(2);
     
    % Linear fitting to get luminance settings (cd/m2).
    x = allCalData{ss}.data.redIntensity;
    nReds = size(allCalData{ss}.dataBC.spdRedBC,2);
    for rr = 1:nReds
        y(rr) = SPDToLuminance(allCalData{ss}.dataBC.spdRedBC(:,rr));
    end 
    p = polyfit(x, y, 1);
    redLuminance = polyval(p, redInputSettings);
    greenLuminance = SPDToLuminance(allCalData{ss}.dataBC.spdGreenBC);
    
    % Calculate mean and std of the output settings.
    meanRedLuminance = mean(redLuminance);
    stdRedLuminance = std(redLuminance);
    
    % Plot raw data (Red).
    plot(ss*ones(1,length(redLuminance)), redLuminance,'ko');
    
    % Plot the reference (Green).
    plot(ss, greenLuminance,'go','markersize',9,'markeredgecolor','k','markerfacecolor',[0.1 0.7 0.1]);
    
    % Plot mean and its standard deviation (Red).
    plot(ss, meanRedLuminance,'ro','markersize',9,'markeredgecolor','k','markerfacecolor','r');
    errorbar(ss, meanRedLuminance, stdRedLuminance, 'r', 'linewidth', 1);
    
    % Save out some results for plotting it again after normalization.
    redLuminanceStruct{ss} = redLuminance;
    greenLuminanceStruct{ss} = greenLuminance;
    greenLuminanceALL(ss) = greenLuminance;
    
    % Print out the progress.
    fprintf('Progress - subject (%d/%d) \n', ss, nSubjects);
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
ylabel('Luminance (cd/m2)','fontsize',fontSize);
xlim([1.5 str2double(numSubjects{end})])
xticks([str2double(numSubjects{1}) : 1 : str2double(numSubjects{end})]);
title('Iso-luminance determination results','fontsize',fontSize);

meanGreenLuminance = mean(greenLuminanceALL);
subtitle(sprintf('This is normalized results for green to be its mean luminance = %.2f (cd/m2)',meanGreenLuminance));

nSubjects = length(numSubjects);
for ss = 1:nSubjects
    % Normalize it.
    luminanceNormalizationFactor = meanGreenLuminance/greenLuminanceStruct{ss};
    greenLuminance = greenLuminanceStruct{ss}*luminanceNormalizationFactor;
    redLuminance = redLuminanceStruct{ss}.*luminanceNormalizationFactor;
    
    % Calculate mean and std of the output settings.
    meanRedLuminance = mean(redLuminance);
    stdRedLuminance = std(redLuminance);
    
    % Plot raw data (Red).
    plot(ss*ones(1,length(redLuminance)), redLuminance,'ko');
    
    % Plot the reference (Green).
    plot(ss, greenLuminance,'go','markersize',9,'markeredgecolor','k','markerfacecolor',[0.1 0.7 0.1]);
    
    % Plot mean and its standard deviation (Red).
    plot(ss, meanRedLuminance,'ro','markersize',9,'markeredgecolor','k','markerfacecolor','r');
    errorbar(ss, meanRedLuminance, stdRedLuminance, 'r', 'linewidth', 1);
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
