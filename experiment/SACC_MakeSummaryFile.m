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
%    09/11/23   smo    - Added sanity check for AUC file.
%    10/12/23   smo    - Added an option to run using either raw data or
%                        processed data without bad contrasts.

%% Initialize.
clear; close all;

%% Set some variables which data to use.
%
% Set this to true if you want to work with the data that we refitted with
% bad points omitted.
FITPFONLYGOODTESTCONTRASTS = true;

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

%% Save out the table.
%
% We will create two separate excel files, one contaning the CS results and
% the other containing the AUC results.
%
% If we use the original results with bad points included, we will do some
% additional analysis and add one more column indicating if the session
% includes bad contrast. If not, we will skip that part and just save out
% the table and stop the routine.
if (SAVETHERESULTS)
    % When using only good data points.
    if (FITPFONLYGOODTESTCONTRASTS)
        range = 'B2';
        testFilenameCS = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_CS.xlsx');
        testFilenameAUC = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
        
        % Make CS summary file.
        writetable(tableCS, testFilenameCS,'Range',range);
        disp('Summary file has been saved successfully! - (CS)');
        
        % Make AUC summary file.
        writetable(tableAUC,testFilenameAUC,'Range',range);
        disp('Summary file has been saved successfully! - (AUC)');
        
        % When using the original dataset with bad points included. It will
        % do additional analysis to sort out the bad sessions.
    elseif ~FITPFONLYGOODTESTCONTRASTS
        
        %% Sanity check for the CSF data.
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
        
        logSensitivitySanityNormal = log10(1/imageContrastNormal);
        logSensitivitySanityHigh = log10(1/imageContrastHigh);
        
        % Get the subject info.
        subjectOptions = unique(tableAUC.Subject);
        nSubjects = length(subjectOptions);
        
        % Read out the test image profile to figure out which test image set was
        % used for 18 cpd per each subject.
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
        testFilename = GetMostRecentFileName(testFiledir,'TestImageProfile.xlsx');
        tableImage = readtable(testFilename);
        
        % We will delete the data in the imageData who didn't complete the study.
        % The subject options who finished is available in 'subjectOptions'.
        subjectToDelete = setdiff(tableImage.Subject,subjectOptions);
        rowToDelete = strcmp(tableImage.Subject,subjectToDelete);
        tableImage(rowToDelete,:) = [];
        
        % Sort the table in the ascending order of spatial frequency.
        tableImage = sortrows(tableImage,'SpatialFrequency','ascend');
        tableImage = sortrows(tableImage,'Subject','ascend');
        
        % Find the subjects who used the high test image set in the experiment.
        testImageMaxContrastHighImageSet = 0.1;
        numSubjectsWithHighImageSet = tableImage.Subject(find(tableImage.TestImageContrastMax == testImageMaxContrastHighImageSet));
        
        % Make a loop to check the sanity per each subject.
        nSessions = size(tableAUC,1);
        for ss = 1:nSessions
            % Get the current number of the subject.
            numSubjectTemp = table2array(tableAUC(ss,'Subject'));
            % Set the criteria differetly if high image set was used.
            if ismember(numSubjectTemp,numSubjectsWithHighImageSet)
                % High image.
                logSensitivitySanity18cpd = logSensitivitySanityHigh;
            else
                % Normal image.
                logSensitivitySanity18cpd = logSensitivitySanityNormal;
            end
            
            % Check if any sensitivity point on CSF is outside of the sanity range.
            tableVarNames = {'LogSensitivity_3cpd','LogSensitivity_6cpd','LogSensitivity_9cpd','LogSensitivity_12cpd','LogSensitivity_18cpd'};
            tableLogSensitivityTemp = table2array(tableAUC(ss,tableVarNames));
            if or(any(tableLogSensitivityTemp(1:4) < logSensitivitySanityNormal),...
                    any(tableLogSensitivityTemp(end) < logSensitivitySanity18cpd))
                idxSessionBad_CSF(ss,1) = 1;
            else
                idxSessionBad_CSF(ss,1) = 0;
            end
        end
        
        % Add a new column to the table saving the findings.
        tableAUC = addvars(tableAUC,idxSessionBad_CSF,'NewVariableNames','IdxSessionBad');
        
        %% Sanity check for the maximum test image contrast used in the experiment.
        nSessions = size(tableImage,1);
        for ss = 1:nSessions
            % Get the current number of the subject.
            numSubjectTemp = table2array(tableImage(ss,'Subject'));
            spatialFrequencyTemp = table2array(tableImage(ss,'SpatialFrequency'));
            testTargetImageContrastMaxTemp = table2array(tableImage(ss,'TestImageContrastMax'));
            
            % Set the criteria differetly for 18 cpd session with high image set.
            if and(spatialFrequencyTemp == 18, testTargetImageContrastMaxTemp == testImageMaxContrastHighImageSet)
                % High image for 18 cpd.
                logSensitivitySanity = logSensitivitySanityHigh;
            else
                % Normal image for the other sessions.
                logSensitivitySanity = logSensitivitySanityNormal;
            end
            
            % Check if the maximum test contrast that used in the experiment was
            % outside of the sanity range. The max contrast value is in
            % 'TestContrasts_8' of the image table.
            testImageContrastMaxTemp = table2array(tableImage(ss,'TestContrasts_8'));
            logSensitiivtyMaxTemp = log10(1/testImageContrastMaxTemp);
            if  (logSensitiivtyMaxTemp < logSensitivitySanity)
                idxSessionBad_maxTestContrast(ss,1) = 1;
            else
                idxSessionBad_maxTestContrast(ss,1) = 0;
            end
        end
        
        % Add the finding to the table.
        tableImage = addvars(tableImage,idxSessionBad_maxTestContrast,'NewVariableNames','IdxSessionBad');
        
        % Get the subjects info of bad sessions.
        subjectSessionBad_maxTestContrast = unique(tableImage.Subject(find(tableImage.IdxSessionBad==1)));
        
        %% Visualize the results what we found from the above.
        figure; clf; hold on;
        spatialFrequencyOptions = [3 6 9 12 18];
        logSpatialFrequencyOptions = log10(spatialFrequencyOptions);
        spatialFrequencyOptionsStr = {'3' '6' '9' '12' '18'};
        
        % Plot the sanity reference.
        plot(log10([1 100]), [logSensitivitySanityNormal logSensitivitySanityNormal],':','linewidth',4,'color',[0 1 0 0.3]);
        plot(log10([1 100]), [logSensitivitySanityHigh logSensitivitySanityHigh],':','linewidth',4,'color',[1 0 0 0.3]);
        
        for ss = 1:nSessions
            % We will use different marker color for Normal image and High image set.
            
            % Update subject name.
            numSubjectTemp = table2array(tableAUC(ss,'Subject'));
            
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
            logSensitivityTemp = table2array(tableAUC(ss,tableVarNames));
            
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
        maxTestContrastsNormal_3cpd = tableImage.TestContrasts_8(find(tableImage.SpatialFrequency == 3));
        maxTestContrastsNormal_6cpd = tableImage.TestContrasts_8(find(tableImage.SpatialFrequency == 6));
        maxTestContrastsNormal_9cpd = tableImage.TestContrasts_8(find(tableImage.SpatialFrequency == 9));
        maxTestContrastsNormal_12cpd = tableImage.TestContrasts_8(find(tableImage.SpatialFrequency == 12));
        
        % 18 cpd. We separate the arrays between normal and high image set.
        idxMaxTestContrasts_18cpd = find(tableImage.SpatialFrequency == 18);
        maxTestContrasts_18cpd = tableImage(idxMaxTestContrasts_18cpd,{'TestImageContrastMax','TestContrasts_8'});
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
        
        % Plot the sanity reference.
        plot(log10([1 100]), [logSensitivitySanityNormal logSensitivitySanityNormal],':','linewidth',4,'color',[0 1 0 0.3]);
        plot(log10([1 100]), [logSensitivitySanityHigh logSensitivitySanityHigh],':','linewidth',4,'color',[1 0 0 0.3]);
        
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
        
        % Print out the info that is outside the sanity range.
        %
        % CSF.
        tableAUC_Bad = tableAUC(find(tableAUC.IdxSessionBad ==1),{'Subject','Filter','LogSensitivity_12cpd','LogSensitivity_18cpd'});
        
        % Find which spatial frequency session was affected.
        for ff = 1:size(tableAUC_Bad,1)
            % Checking 12 cpd.
            if (tableAUC_Bad.LogSensitivity_12cpd(ff) < logSensitivitySanityNormal)
                whichSF_tableAUC_Bad(ff,1) = 12;
                % If not 12 cpd, it will be 18 cpd, so set it as 18 cpd.
            else
                whichSF_tableAUC_Bad(ff,1) = 18;
            end
        end
        
        % Will add one more column to the table.
        tableAUC_Bad = addvars(tableAUC_Bad,whichSF_tableAUC_Bad,'NewVariableNames','Spatial Frequency (Bad)');
        
        % Max image contrast.
        tableImage_Bad = tableImage(find(tableImage.IdxSessionBad == 1),{'Subject','SpatialFrequency','TestContrasts_8'});
        tableImage_Bad = sortrows(tableImage_Bad,'TestContrasts_8');
        tableImage_Bad = sortrows(tableImage_Bad,'SpatialFrequency');
        
        %% Finally, here we save out the updated table.
        %
        % We will create two separate excel files, one contaning the CS results and
        % the other containing the AUC results.
        range = 'B2';
        testFilenameCS = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_CS.xlsx');
        testFilenameAUC = fullfile(testFiledir, 'SACC_Experiment_Results_Summary_AUC.xlsx');
        
        % Make CS summary file.
        writetable(tableCS, testFilenameCS,'Range',range);
        disp('Summary file has been saved successfully! - (CS)');
        
        % Make AUC summary file.
        writetable(tableAUC,testFilenameAUC,'Range',range);
        disp('Summary file has been saved successfully! - (AUC)');
    end
end
