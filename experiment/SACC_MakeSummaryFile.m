%% SACC_MakeSummaryFile
%
% This creates a text summary file by merging CS and AUC calculation
% results for all subjects. The created file should basically contain all
% results acquired from this project.
%
% See also:
%    t_CSFGeneratorAnalyze, SACC_FitCSF

% History:
%    03/10/23   smo    - Started on it.
%    04/24/23   smo    - Now making CS and AUD files separately.
%    09/11/23   smo    - Added checking good contrasts for AUC file.
%    10/12/23   smo    - Added an option to run using either raw data or
%                        processed data without bad contrasts.

%% Initialize.
clear; close all;

%% Set some variables which data to use.
%
% Set this to true if you want to work with the data that we refitted with
% bad points omitted.
FITPFONLYGOODTESTCONTRASTS = false;

% Set directory differently to save.
if (FITPFONLYGOODTESTCONTRASTS)
    whichPref = 'SACCAnalysisFinal';
else
    whichPref = 'SACCAnalysis';
end

% Decide if you want to save out the table made from this routine.
SAVETHERESULTS = false;

%% Get available subjects.
if ispref('SpatioSpectralStimulator',whichPref)
    testFiledir = fullfile(getpref('SpatioSpectralStimulator',whichPref));
end
testFileList = dir(testFiledir);

% Subject name.
for tt = 1:length(testFileList)
    testFilenameList{tt}  = testFileList(tt).name;
end
idxSubjectName = find(str2double(testFilenameList)>0);
subjectNameOptions = testFilenameList(idxSubjectName);

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
t_CS = vertcat(dataCS{:});

% Update the numbering.
nRows = size(t_CS,1);
t_CS{:,'No'} = linspace(1,nRows,nRows)';

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
t_AUC = vertcat(dataAUC{:});

% Update the numbering.
nRows = size(t_AUC,1);
t_AUC{:,'No'} = linspace(1,nRows,nRows)';

disp('AUC data has been merged successfully!');

%% Save out the table.
%
% We will create two separate excel files, one contaning the CS results and
% the other containing the AUC results.
%
% If we use the original results with bad points included (when set
% 'FITPFONLYGOODTESTCONTRASTS' to false), we will do some additional
% analysis and add one more column indicating if the session includes bad
% contrast. If not, we will skip that part and just save out the table and
% stop the routine.
%
% When using only good data points. We will save the table and stop the
% routine.
if (FITPFONLYGOODTESTCONTRASTS)
    if (SAVETHERESULTS)
        range = 'B2';
        testFilenameCS = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_CS.xlsx');
        testFilenameAUC = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
        
        % Make CS summary file.
        writetable(t_CS, testFilenameCS,'Range',range);
        disp('Summary file has been saved successfully! - (CS)');
        
        % Make AUC summary file.
        writetable(t_AUC,testFilenameAUC,'Range',range);
        disp('Summary file has been saved successfully! - (AUC)');
    end
    
    % When using the original dataset with bad points included. It will
    % do additional analysis to sort out the bad sessions.
elseif ~FITPFONLYGOODTESTCONTRASTS
    %% Check if the CSF data is within a good contrast range.
    %
    % This will be run only if we use the original data that used all the test
    % contrasts including the bad ones (as of 10/12/23).
    %
    % As we found out that we used the test images using the Standard method,
    % here we want to check if the CSF data was fitted within the 'good'
    % contrast range of the test images.
    %
    % These cone contrasts are the marginal contrasts that can generate the
    % desired contrast. If a test image contrast had a higher contrast than
    % this, its contrast would be not perfectly modulated as desired.
    %
    % These values were found using the validation code
    % SpectralCalCheck_ver2.m. This routine searches the affected sessions
    % based on these information - coneContrastNormal and coneContrastHigh - so
    % we can update these as we want.
    coneContrastNormal = [-0.0427 0.0369 -0.0028];
    coneContrastHigh = [-0.0579 0.0590 -0.0003];
    
    imageContrastNormal = sqrt(sum(coneContrastNormal.^2));
    imageContrastHigh = sqrt(sum(coneContrastHigh.^2));
    
    marginalLogSensitivityNormal = log10(1/imageContrastNormal);
    marginalLogSensitivityHigh = log10(1/imageContrastHigh);
    
    % Get the subject info.
    subjectOptions = unique(t_AUC.Subject);
    nSubjects = length(subjectOptions);
    
    % Read out the test image profile to figure out which test image set was
    % used at 18 cpd session per each subject.
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
    testFilename = GetMostRecentFileName(testFiledir,'TestImageProfile.xlsx');
    t_Image = readtable(testFilename);
    
    % We will delete the data in the imageData who didn't complete the study.
    % The subject options who finished is available in 'subjectOptions'.
    subjectToDelete = setdiff(t_Image.Subject,subjectOptions);
    rowToDelete = strcmp(t_Image.Subject,subjectToDelete);
    t_Image(rowToDelete,:) = [];
    
    % Sort the table in the ascending order of spatial frequency.
    t_Image = sortrows(t_Image,'SpatialFrequency','ascend');
    t_Image = sortrows(t_Image,'Subject','ascend');
    
    % Find the subjects who used the high test image set in the experiment.
    maxContrastHighImageSet = 0.1;
    numSubjectsWithHighImageSet = t_Image.Subject(find(t_Image.TestImageContrastMax == maxContrastHighImageSet));
    
    % Make a loop to check if test contrasts are in a good range per each
    % subject.
    nSessions = size(t_AUC,1);
    for ss = 1:nSessions
        % Get the current number of the subject.
        numSubjectTemp = table2array(t_AUC(ss,'Subject'));
        % Set the criteria differetly if high image set was used.
        if ismember(numSubjectTemp,numSubjectsWithHighImageSet)
            % High image.
            marginalLogSensitivity18cpd = marginalLogSensitivityHigh;
        else
            % Normal image.
            marginalLogSensitivity18cpd = marginalLogSensitivityNormal;
        end
        
        % Get the sensitivity data per subject.
        tableVarNames = {'LogSensitivity_3cpd','LogSensitivity_6cpd','LogSensitivity_9cpd','LogSensitivity_12cpd','LogSensitivity_18cpd'};
        logSensitivityPerSessionTemp = table2array(t_AUC(ss,tableVarNames));
        
        % Check if any sensitivity point on the CSF is outside of the
        % good range.
        %
        % For some subjects at 18 cpd session, we used high contrast
        % image, so we set the marginal contrast (sensitivity)
        % differently.
        if or(any(logSensitivityPerSessionTemp(1:4) < marginalLogSensitivityNormal),...
                any(logSensitivityPerSessionTemp(end) < marginalLogSensitivity18cpd))
            idxSessionBad_CSF(ss,1) = 1;
        else
            idxSessionBad_CSF(ss,1) = 0;
        end
    end
    
    % Add a new column to the table saving the findings.
    t_AUC = addvars(t_AUC,idxSessionBad_CSF,'NewVariableNames','IdxSessionBad');
    
    % Print out the info that is outside the good contrast range.
    t_AUCBad = t_AUC(find(t_AUC.IdxSessionBad ==1),{'Subject','Filter','LogSensitivity_12cpd','LogSensitivity_18cpd'});
    
    % Find which spatial frequency session was affected.
    for ff = 1:size(t_AUCBad,1)
        % Checking 12 cpd.
        if (t_AUCBad.LogSensitivity_12cpd(ff) < marginalLogSensitivityNormal)
            whichSF_tableAUC_Bad(ff,1) = 12;
            % If not 12 cpd, it will be 18 cpd, so set it as 18 cpd.
        else
            whichSF_tableAUC_Bad(ff,1) = 18;
        end
    end
    
    % Will add one more column to the table.
    t_AUCBad = addvars(t_AUCBad,whichSF_tableAUC_Bad,'NewVariableNames','Spatial Frequency (Bad)');
    
    %% Check if the maximum test image contrast used in the experiment is in a good range.
    %
    % We will make a loop to check each session. The 'nSessions' should be
    % 160 (32 subjects x 5 spatial frequencies).
    nSessions = size(t_Image,1);
    for ss = 1:nSessions
        % Get the current subject info.
        numSubjectTemp = table2array(t_Image(ss,'Subject'));
        spatialFrequencyTemp = table2array(t_Image(ss,'SpatialFrequency'));
        testTargetImageContrastMaxTemp = table2array(t_Image(ss,'TestImageContrastMax'));
        
        % Set the criteria differetly for 18 cpd session with high image set.
        if and(spatialFrequencyTemp == 18, testTargetImageContrastMaxTemp == maxContrastHighImageSet)
            % High image for 18 cpd.
            margianlLogSensitivity = marginalLogSensitivityHigh;
        else
            % Normal image for the other sessions.
            margianlLogSensitivity = marginalLogSensitivityNormal;
        end
        
        % Check if the maximum test contrast that used in the experiment was
        % outside of the good contrast range. The max contrast value is in
        % 'TestContrasts_8' of the image table.
        testImageContrastMaxTemp = table2array(t_Image(ss,'TestContrasts_8'));
        logSensitiivtyMaxTemp = log10(1/testImageContrastMaxTemp);
        if  (logSensitiivtyMaxTemp < margianlLogSensitivity)
            idxSessionBad_maxTestContrast(ss,1) = 1;
        else
            idxSessionBad_maxTestContrast(ss,1) = 0;
        end
    end
    
    % Add the finding to the table.
    t_Image = addvars(t_Image,idxSessionBad_maxTestContrast,'NewVariableNames','IdxSessionBad');
    
    % Find the sessions with bad contrast presented.
    t_ImageBad = t_Image(find(t_Image.IdxSessionBad == 1),{'Subject','SpatialFrequency','TestContrasts_8'});
    t_ImageBad = sortrows(t_ImageBad,'TestContrasts_8');
    t_ImageBad = sortrows(t_ImageBad,'SpatialFrequency');
    
    %% Visualize the results what we found from the above.
    figure; clf; hold on;
    spatialFrequencyOptions = [3 6 9 12 18];
    logSpatialFrequencyOptions = log10(spatialFrequencyOptions);
    spatialFrequencyOptionsStr = {'3' '6' '9' '12' '18'};
    
    % Plot the good contrast range as a reference.
    plot(log10([1 100]), [marginalLogSensitivityNormal marginalLogSensitivityNormal],':','linewidth',4,'color',[0 1 0 0.3]);
    plot(log10([1 100]), [marginalLogSensitivityHigh marginalLogSensitivityHigh],':','linewidth',4,'color',[1 0 0 0.3]);
    
    for ss = 1:nSessions
        % We will use different marker color for Normal image and High image set.
        
        % Update subject name.
        numSubjectTemp = table2array(t_AUC(ss,'Subject'));
        
        % Set marker color differently over test image set.
        markerColorNormal = 'g';
        markerColorHigh = 'r';
        if ismember(numSubjectTemp,numSubjectsWithHighImageSet)
            markerColor = markerColorHigh;
        else
            markerColor = markerColorNormal;
        end
        
        % Get CSF curve values.
        tableVarNames = {'LogSensitivity_3cpd','LogSensitivity_6cpd','LogSensitivity_9cpd','LogSensitivity_12cpd','LogSensitivity_18cpd'};
        logSensitivityTemp = table2array(t_AUC(ss,tableVarNames));
        
        % Plot it - 3, 6, 9, 12 cpd.
        % We only used a normal image set for this spatial frequencies, so it
        % will be all same colors.
        plot(logSpatialFrequencyOptions(1:4),logSensitivityTemp(1:4),'ko','markerfacecolor',markerColorNormal);
        
        % Plot it - 18 cpd.
        % For 18 cpd, we will plot it in different colors per different image
        % set we used (normal / high).
        plot(logSpatialFrequencyOptions(end),logSensitivityTemp(end),'ko','markerfacecolor',markerColor);
    end
    
    % Set axis and stuffs.
    % This is the same format as our analysis.
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Log Contrast Sensitivity','fontsize',15);
    xticks(logSpatialFrequencyOptions);
    xticklabels(spatialFrequencyOptionsStr);
    xlim([min(log10(spatialFrequencyOptions)) max(log10(spatialFrequencyOptions))]);
    ylim(log10([1 600]));
    yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
    yticks(yaxisRange);
    ytickformat('%.2f');
    title(sprintf('CSF data points (N=%d)',nSessions), 'fontsize',15);
    f = flip(get(gca,'children'));
    legend(f([1 2 3 100]), 'Ref-Normal','Ref-High','CSF (Normal)', 'CSF (High)','location','southeast','fontsize',15);
    
    %% We will plot another figure in the same format as above, but now with the
    % highest contrast image presented in the experiment.
    figure; clf; hold on;
    
    % 3, 6, 9, and 12 cpd (All normal image set).
    maxTestContrastsNormal_3cpd = t_Image.TestContrasts_8(find(t_Image.SpatialFrequency == 3));
    maxTestContrastsNormal_6cpd = t_Image.TestContrasts_8(find(t_Image.SpatialFrequency == 6));
    maxTestContrastsNormal_9cpd = t_Image.TestContrasts_8(find(t_Image.SpatialFrequency == 9));
    maxTestContrastsNormal_12cpd = t_Image.TestContrasts_8(find(t_Image.SpatialFrequency == 12));
    
    % 18 cpd. We separate the arrays between normal and high image set.
    idxMaxTestContrasts_18cpd = find(t_Image.SpatialFrequency == 18);
    maxTestContrasts_18cpd = t_Image(idxMaxTestContrasts_18cpd,{'TestImageContrastMax','TestContrasts_8'});
    maxTestContrastOptions_18cpd = unique(maxTestContrasts_18cpd.TestImageContrastMax);
    
    maxTestContrastNormal = min(maxTestContrastOptions_18cpd);
    maxTestContrastHigh = max(maxTestContrastOptions_18cpd);
    
    maxTestContrastsNormal_18cpd = maxTestContrasts_18cpd.TestContrasts_8(maxTestContrasts_18cpd.TestImageContrastMax == maxTestContrastNormal);
    maxTestContrastsHigh_18cpd = maxTestContrasts_18cpd.TestContrasts_8(maxTestContrasts_18cpd.TestImageContrastMax == maxTestContrastHigh);
    
    % Convert to log sensitivity.
    logSensitivityMaxTestContrastsNormal_3cpd = log10(1./maxTestContrastsNormal_3cpd);
    logSensitivityMaxTestContrastsNormal_6cpd = log10(1./maxTestContrastsNormal_6cpd);
    logSensitivityMaxTestContrastsNormal_9cpd = log10(1./maxTestContrastsNormal_9cpd);
    logSensitivityMaxTestContrastsNormal_12cpd = log10(1./maxTestContrastsNormal_12cpd);
    logSensitivityMaxTestContrastsNormal_18cpd = log10(1./maxTestContrastsNormal_18cpd);
    logSensitivityMaxTestContrastsHigh_18cpd = log10(1./maxTestContrastsHigh_18cpd);
    
    % Plot the good contrast range as a reference.
    plot(log10([1 100]), [marginalLogSensitivityNormal marginalLogSensitivityNormal],':','linewidth',4,'color',[0 1 0 0.3]);
    plot(log10([1 100]), [marginalLogSensitivityHigh marginalLogSensitivityHigh],':','linewidth',4,'color',[1 0 0 0.3]);
    
    % Plot the sensitivity of the maximum test image contrast.
    %
    % Normal image.
    plot(logSpatialFrequencyOptions(1),logSensitivityMaxTestContrastsNormal_3cpd,'ko','markerfacecolor',markerColorNormal);
    plot(logSpatialFrequencyOptions(2),logSensitivityMaxTestContrastsNormal_6cpd,'ko','markerfacecolor',markerColorNormal);
    plot(logSpatialFrequencyOptions(3),logSensitivityMaxTestContrastsNormal_9cpd,'ko','markerfacecolor',markerColorNormal);
    plot(logSpatialFrequencyOptions(4),logSensitivityMaxTestContrastsNormal_12cpd,'ko','markerfacecolor',markerColorNormal);
    plot(logSpatialFrequencyOptions(5),logSensitivityMaxTestContrastsNormal_18cpd,'ko','markerfacecolor',markerColorNormal);
    
    % High image.
    plot(logSpatialFrequencyOptions(5),logSensitivityMaxTestContrastsHigh_18cpd,'ko','markerfacecolor',markerColorHigh);
    
    % Set axis and stuffs.
    % This is the same format as our analysis.
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Log Contrast Sensitivity','fontsize',15);
    xticks(logSpatialFrequencyOptions);
    xticklabels(spatialFrequencyOptionsStr);
    xlim([min(log10(spatialFrequencyOptions)) max(log10(spatialFrequencyOptions))]);
    ylim(log10([1 600]));
    yaxisRange = log10([10, 100, 200, 300, 400, 500, 600]);
    yticks(yaxisRange);
    ytickformat('%.2f');
    title(sprintf('Max test image contrast (N=%d)',nSessions), 'fontsize',15);
    f = flip(get(gca,'children'));
    legend(f([1 2 3 end]), 'Ref-Normal','Ref-High','CSF (Normal)', 'CSF (High)','location','southeast','fontsize',15);
      
    %% Finally, here we save out the updated table.
    %
    % We will create two separate excel files, one contaning the CS results and
    % the other containing the AUC results.
    if (SAVETHERESULTS)
        range = 'B2';
        testFilenameCS = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_CS.xlsx');
        testFilenameAUC = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
        
        % Make CS summary file.
        writetable(t_CS, testFilenameCS,'Range',range);
        disp('Summary file has been saved successfully! - (CS)');
        
        % Make AUC summary file.
        writetable(t_AUC,testFilenameAUC,'Range',range);
        disp('Summary file has been saved successfully! - (AUC)');
    end
end
