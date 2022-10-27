% t_CSFGeneratorAnalyze.
%
% Syntax:
%    t_CSFGeneratorAnalyze
%
% Description:
%    This is for fitting PF to the data acquired from the experiments.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    t_CSFGeneratorExperiment.

% History:
%    02/28/22  smo          - Started on it.
%    03/14/22  smo          - Added a plotting option to make different
%                             marker size over different number of trials.
%    03/16/22  smo          - Added an option to check if Adaptive mode
%                             works fine.
%    06/26/22  smo          - Added a part plotting CSF curves.
%    10/03/22  smo          - Now fitting all spatial frequency at once and
%                             give CSF curve at the same time.
%    10/26/22  smo          - Added an option to fit all available data at
%                             once.

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
CHECKADAPTIVEMODE = false;
PF = 'weibull';
paramsFree = [1 1 0 1];

olderDate = 0;
SUBPLOT = true;
axisLog = true;
addInitialThresholdEst = true;
addQuestFit = true;
addLegend = false;

%% Find available data here.
%
% You can either fit all available data at once or choose one subject's
% data to fit.
FITALLATONCE = true;

if (FITALLATONCE)
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        testFileList = dir(testFiledir);
        
        % Find available data of subject and spatial frequency.
        %
        % Subject name.
        for tt = 1:length(testFileList)
            testFilenameList{tt}  = testFileList(tt).name;
        end
        idxSubjectName = find(str2double(testFilenameList)>0);
        subjectNameOptions = testFilenameList(idxSubjectName);
        
        % Spatial frequency.
        nSubjectName = length(subjectNameOptions);
        for ss = 1:nSubjectName
            dirTemp = fullfile(testFiledir,subjectNameOptions{ss});
            dataList = dir(dirTemp);
            
            numStartData = 3;
            for dd = 1:length(dataList)-numStartData+1
                idxData = numStartData+dd-1;
                spatialFrequencyOptions{dd,ss} = dataList(idxData).name;
            end
        end
    end
else
    % Load single subject to fit one by one.
    subjectNameOptions = {'012'};
    spatialFrequencyOptions = {'3'};
end

%% Show the progress of the experiment.
SHOWPROGRESS = true;

if (SHOWPROGRESS)
    figure; clf;
    
    % Set background axes. This is just for texts.
    main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
    
    % Put visit number text on the top.
    nVisits = 5;
    barPositionHorzStart = 0.15;
    barPositionHorzEnd = 0.8;
    barPositionVert = 0.8;
    for vv = 1:nVisits
        visitOptions = linspace(3,7,nVisits);
        visitStr = sprintf('Visit %d', visitOptions(vv));
        textIntervalHorzOptions = linspace(barPositionHorzStart+0.05,barPositionHorzEnd,nVisits);
        text(textIntervalHorzOptions(vv), barPositionVert+0.05, visitStr, 'Parent', main);
    end
    
    % Use a red annotation rectangle as background, and overlay a
    % green annotation rectangle on top.
    nSubjects = length(subjectNameOptions);
    for ss = 1:nSubjects
        % Subejct.
        textLocationVertOptions = sort(linspace(0.2,0.8,nSubjects),'descend');
        textLocationVert = textLocationVertOptions(ss);
        text(0.03, textLocationVert, sprintf('Subject %s',subjectNameOptions{ss}), 'Parent', main)
        
        % Fill the bar based on the completed number of spatial frequency
        % (so, number of visits) per each subject.
        barWidth = 0.02;

        numVisitsCompleted = spatialFrequencyOptions(:,ss);
        numVisitsCompleted = numVisitsCompleted(find(~cellfun(@isempty,numVisitsCompleted)));
        numLevelProgress = length(numVisitsCompleted);
        
        bgBarColor = annotation('rectangle', ...
            [barPositionHorzStart textLocationVert barPositionHorzEnd barWidth],...
            'EdgeColor','black', 'FaceColor', 'white');
        frontBarColor = annotation('rectangle', ...
            [barPositionHorzStart textLocationVert numLevelProgress*0.15 barWidth],...
            'EdgeColor','None', 'FaceColor', 'blue');
    end
    
    % Save the progress plot.
    SAVEPROGRESSPLOT = true;
    if (SAVEPROGRESSPLOT)
        
    end
end



%% Load data and PF fitting.
nSubjects = length(subjectNameOptions);

for ss = 1:nSubjects
    % Set target subject.
    subjectName = subjectNameOptions{ss};
    
    % Set available spatial frequency data for the subject.
    sineFreqCyclesPerDeg = spatialFrequencyOptions(:,ss);
    % Here we remove empty cells. It shows a cell empty when a subject does
    % not have all possible spatial frequency data.
    sineFreqCyclesPerDeg = sineFreqCyclesPerDeg(...
        find(~cellfun(@isempty,sineFreqCyclesPerDeg)));
    
    nSineFreqCyclesPerDeg = length(sineFreqCyclesPerDeg);
    
    % Set all possible filters.
    filterOptions = {'A', 'B', 'C', 'D', 'E'};
    nFilters = length(filterOptions);
    
    % Set the size of the subplot. Each figure will contain all filters
    % data, so there will be a total of five subplots in one figure.
    sizeSubplot = [2 3];
    
    for dd = 1:nSineFreqCyclesPerDeg
        
        % As spatial frequency info is stored in string form, we extract
        % and convert it to double.
        sineFreqCyclesPerDegStr = sineFreqCyclesPerDeg{dd};
        sineFreqCyclesPerDegTemp = sscanf(sineFreqCyclesPerDegStr,'%d');
        
        figure; clf;
        set(gcf,'position',[0,0,1920,1080]);
        
        for ff = 1:nFilters
            % Set target filter.
            whichFilter = filterOptions{ff};
            
            % Load the experiment data.
            if (ispref('SpatioSpectralStimulator','SACCData'))
                testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),...
                    subjectName,append(num2str(sineFreqCyclesPerDegTemp),'_cpd'));
                testFilename = GetMostRecentFileName(testFiledir,...
                    sprintf('CS_%s_%d_cpd_%s',subjectName,sineFreqCyclesPerDegTemp,whichFilter),'olderDate',olderDate);
                theData = load(testFilename);
            else
                error('Cannot find data file');
            end
            
            nDataContrastRange = 1;
            for cc = 1:nDataContrastRange
                % Load the contrast range data.
                testFilename = GetMostRecentFileName(testFiledir,...
                    sprintf('ContrastRange_%s_%d',subjectName,sineFreqCyclesPerDegTemp), 'olderDate',olderDate+cc-1);
                theContrastData = load(testFilename);
                
                % Extract the threshold from the initial measurements.
                thresholdInitial(cc,:) = theContrastData.preExpDataStruct.thresholdFoundRawLinear;
            end
            
            % Pull out the data here.
            nTrials = theData.estimator.nRepeat;
            [stimVec, responseVec, structVec] = combineData(theData.estimator);
            
            % Here we will plot the PF fitting graph.
            %
            % Set marker size here. We will use this number to plot the results to
            % have different marker size according to the number of trials. Here
            % we used the same method to decide the size of each marker as
            % 'thresholdMLE' does.
            stimVal = unique(stimVec);
            pCorrect = zeros(1,length(stimVal));
            for idx = 1:length(stimVal)
                prop = responseVec(stimVec == stimVal(idx));
                pCorrect(idx) = sum(prop) / length(prop);
                pointSize(idx) = 10 * 100 / length(stimVec) * length(prop);
            end
            
            thresholdCriterion = 0.81606;
            [threshold, para, dataOut] = theData.estimator.thresholdMLE(...
                'thresholdCriterion', thresholdCriterion, 'returnData', true);
            thresholdsQuest(ff) = threshold;
            
            % Set the contrast levels in linear unit.
            examinedContrastsLinear = 10.^dataOut.examinedContrasts;
            
            % Set if you want to add Quest fit in the results.
            if (addQuestFit)
                questPara = para;
            else
                questPara = [];
            end
            
            % PF fitting here.
            if (SUBPLOT)
                subplot(sizeSubplot(1),sizeSubplot(2),ff); hold on;
            end
            [paramsFitted(:,ff)] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
                'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
                'newFigureWindow', ~SUBPLOT, 'pointSize', pointSize, 'axisLog', axisLog,...
                'questPara', questPara,'addLegend',false);
            subtitle(sprintf('%d cpd / Filter = %s',sineFreqCyclesPerDegTemp,whichFilter),'fontsize', 15);
            
            % Add initial threhold to the plot.
            if (addInitialThresholdEst)
                for cc = 1:nDataContrastRange
                    % Plot it on log space if you want.
                    if(axisLog)
                        thresholdInitial(cc,:) = log10(thresholdInitial(cc,:));
                    end
                    
                    % Plot it here.
                    plot([thresholdInitial(cc,1) thresholdInitial(cc,1)], [0 1], 'b-', 'linewidth',3);
                    plot([thresholdInitial(cc,2) thresholdInitial(cc,2)], [0 1], 'g--', 'linewidth',3);
                    plot([thresholdInitial(cc,3) thresholdInitial(cc,3)], [0 1], 'b-', 'linewidth',3);
                    plot([thresholdInitial(cc,4) thresholdInitial(cc,4)], [0 1], 'g--', 'linewidth',3);
                end
            end
            
            if (addQuestFit)
                legend('Data','PF-fit','PF-Threshold','Quest-fit','ThresholdEst from high','ThresholdEst from low',...
                    'FontSize', 9, 'location', 'southeast');
            else
                legend('Data','PF-fit','PF-Threshold','ThresholdEst from high','ThresholdEst from low',...
                    'FontSize', 9, 'location', 'southeast');
            end
            
            % Set the range for the x-axis.
            xlim([-3.3 -1]);
            
            % Clear the pointsize for next plot.
            clear pointSize;
        end
        
        % Print out the progress.
        fprintf('\t Fitting progress - Subject (%d/%d) / Spatial frequency (%d/%d) \n',ss, nSubjects, dd, nSineFreqCyclesPerDeg);
        
        %% Save the results.
        SAVETHEPLOT = true;
        
        if (SAVETHEPLOT)
            if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),...
                    subjectName,append(num2str(sineFreqCyclesPerDegTemp),'_cpd'));
                
                % Make folder with subject name if it does not exist.
                if ~exist(testFiledir, 'dir')
                    mkdir(testFiledir);
                end
                
                % Save the plot.
                testFilename = fullfile(testFiledir,...
                    sprintf('CS_%s_%d_cpd',subjectName,sineFreqCyclesPerDegTemp));
                testFileFormat = '.tiff';
                saveas(gcf,append(testFilename,testFileFormat));
                fprintf('\t Plot has been saved successfully! \n');
            end
        end
    end
end

%% Here we check if the Adaptive mode works fine if you want.
if (CHECKADAPTIVEMODE)
    figure; clf; hold on;
    for ss = 1:nSineFreqCyclesPerDeg
        
        % Set target spatial frequency.
        sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
        
        % Load the data.
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),...
                subjectName,append(num2str(sineFreqCyclesPerDegTemp),'_cpd'));
            testFilename = GetMostRecentFileName(testFiledir,...
                'CS','olderDate',olderDate);
            theData = load(testFilename);
        else
            error('Cannot find data file');
        end
        
        % Set variables here.
        nTrial = theData.estimator.nTrial;
        testTrials = linspace(1,nTrial,nTrial);
        testContrasts = 10.^stimVec;
        testPerformances = {structVec.outcome};
        
        % Set the color of the data point according to the subject's
        % performance. We will differenciate the marker point based on
        % each response either correct of incorrect.
        %
        % In testPerformances, correct response is allocated to 2 and
        % incorrect to 1. We will differenciate the marekr color based on
        % these.
        responseCorrect   = 2;
        responseIncorrect = 1;
        markerColorCorrect   = 5;
        markerColorIncorrect = 10;
        
        for tt = 1:nTrial
            markerFaceColor(1,tt) = testPerformances{tt};
            markerFaceColor(find(markerFaceColor == responseCorrect)) = markerColorCorrect;
            markerFaceColor(find(markerFaceColor == responseIncorrect)) = markerColorIncorrect;
        end
        
        if (SUBPLOT)
            subplot(sizeSubplot(1),sizeSubplot(2),ss); hold on;
        end
        
        % Plot it.
        markerSize = 40;
        scatter(testTrials, testContrasts, markerSize, markerFaceColor, 'filled', 'MarkerEdgeColor', zeros(1,3));
        xlabel('Number of Trials', 'fontsize', 15);
        ylabel('Contrast', 'fontsize', 15);
        title(sprintf('%d cpd',sineFreqCyclesPerDegTemp),'fontsize', 15);
        legend('Blue = correct / Yellow = incorrect','location','northeast','fontsize',15);
        
        clear markerFaceColor;
    end
end

%% Plot the CSF curve here.
CSFCURVE = false;
SUBPLOTCSF = false;

if (CSFCURVE)
    % Export the threshold data.
    thresholds = paramsFitted(1,:);
    
    % Convert NaN to 0 here.
    for tt = 1:length(thresholds)
        if isnan(thresholds(tt))
            thresholds(tt) = 0;
        end
    end
    
    % Calculate sensitivity.
    sensitivityLinear = 1./thresholds;
    sensitivityLog = log10(sensitivityLinear);
    sineFreqCyclesPerDegLog = log10(sineFreqCyclesPerDeg);
    
    sensitivityQuestLinear = 1./(10.^thresholdsQuest);
    sensitivityQuestLog = log10(sensitivityQuestLinear);
    
    % Decide the plot on either subplot or separate figure.
    if (SUBPLOTCSF)
        subplot(sizeSubplot(1), sizeSubplot(2), 6); hold on;
    else
        figure; clf; hold on;
    end
    
    % Plot PF fit CSF curve.
    plot(sineFreqCyclesPerDegLog, sensitivityLog, 'r.-','markersize',20,'linewidth',2);
    % Add Quest fit CSF curve.
    if (addQuestFit)
        plot(sineFreqCyclesPerDegLog, sensitivityQuestLog, 'k.--','markersize',20,'linewidth',2);
    end
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast Sensitivity','fontsize',15);
    xticks(sineFreqCyclesPerDegLog);
    xticklabels(sineFreqCyclesPerDeg);
    yticks(sort(sensitivityLog));
    yticklabels(sort(round(sensitivityLinear)));
    title('CSF curve','fontsize',15);
    % Add legend.
    if (addQuestFit)
        legend(append(subjectName,'-PF'),append(subjectName,'-Quest'),'fontsize',15);
    else
        legend(append(subjectName,'-PF'),'fontsize',15);
    end
end
