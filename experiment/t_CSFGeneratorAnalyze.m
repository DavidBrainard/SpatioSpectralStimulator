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
%    11/14/22  dhb          - Bootstrapping
%    11/17/22  dhb          - Two pass fitting.  Needs a code cleaning pass.

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
CHECKADAPTIVEMODE = false;

FITALLATONCE = true;
SAVETHEPLOT = true;
RECORDTESTIMAGEPROFILE = true;

% Psychometric function information.  
% The paramsFree slope field is overridden
% when we pass a range of slopes explicitly.
PF = 'weibull';
paramsFree = [1 1 0 1];

% Allows this to analyze older data files.  0
% means most recent, 1 next most recent, etc.
% For actual (non-pilot) subjects, there should
% be only 1 data file in each top level data directory
% and this should be 0.
olderDate = 0;

% Control things about the plots
SUBPLOT = true;
axisLog = true;
addInitialThresholdEst = true;
addQuestFit = false;

% Slope (aka beta) values.  Fits are done on this
% discrete set of slopes, with best slope chosen.
% Set to empty to allow the slope to go free, according
% to paramsFree vector above.
%
% An initial fit is done on slopeValList, and then a second
% pass for each session (subject/sf) with slopes limited
% by what was found in the initial fit.  The parameter
% slopeRangeLogUnits affects the expansion of the determined
% range for the second fit.
minSlope = 0.1;
maxSlope = 10;
nSlopes = 20;
slopeValList = 10.^linspace(log10(minSlope),log10(maxSlope),nSlopes);
slopeRangeLogUnits = 0.2;

% Bootstrap info
BOOTSTRAP_RAWFITS = true;
nBootstraps = 100;
bootConfInterval = 0.8;

% Threshold criterion
thresholdCriterion = 0.81606;

%% Find available data here.
%
% You can either fit all available data at once or choose one subject's
% data to fit.
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
    subjectNameOptions = {'003'};
    spatialFrequencyOptions = {'6'};
end

%% Show the progress of the experiment.
if (FITALLATONCE)
    SHOWPROGRESS = true;
    if (SHOWPROGRESS)
        figure; clf;
        
        % Set background axes. This is just for texts.
        main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
        
        % Put visit number text on the top.
        numVisitStart = 3;
        numVisitEnd = 7;
        nVisits = numVisitEnd-numVisitStart+1;
        barPositionHorzStart = 0.15;
        barPositionHorzEnd = 0.8;
        barPositionVert = 0.8;
        for vv = 1:nVisits
            visitOptions = linspace(numVisitStart,numVisitEnd,nVisits);
            visitStr = sprintf('Visit %d', visitOptions(vv));
            textIntervalHorzOptions = linspace(barPositionHorzStart+0.05,barPositionHorzEnd+0.05,nVisits);
            text(textIntervalHorzOptions(vv), barPositionVert+0.05, visitStr, 'Parent', main);
        end
        
        % Use a red annotation rectangle as background, and overlay a
        % green annotation rectangle on top.
        nSubjectsProgress = length(subjectNameOptions);
        for ss = 1:nSubjectsProgress
            % Subejct.
            textLocationVertOptions = sort(linspace(0.2,0.8,nSubjectsProgress),'descend');
            textLocationVert = textLocationVertOptions(ss);
            text(0.03, textLocationVert, sprintf('Subject %s',subjectNameOptions{ss}), 'Parent', main)
            
            % Fill the bar based on the completed number of spatial frequency
            % (so, number of visits) per each subject.
            barWidth = 0.02;
            
            numVisitsCompleted = spatialFrequencyOptions(:,ss);
            numVisitsCompleted = numVisitsCompleted(find(~cellfun(@isempty,numVisitsCompleted)));
            numLevelProgress = length(numVisitsCompleted);
            barWidthHorzOneLevelProgress = 0.16;
            
            bgBarColor = annotation('rectangle', ...
                [barPositionHorzStart textLocationVert barPositionHorzEnd barWidth],...
                'EdgeColor','black', 'FaceColor', 'white');
            frontBarColor = annotation('rectangle', ...
                [barPositionHorzStart textLocationVert numLevelProgress*barWidthHorzOneLevelProgress barWidth],...
                'EdgeColor','None', 'FaceColor', 'blue');
        end
        
        % Save the progress plot.
        SAVEPROGRESSPLOT = true;
        if (SAVEPROGRESSPLOT)
            if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
                testFilename = fullfile(testFiledir,'Experiment_Progress');
                testFileFormat = '.tiff';
                saveas(gcf,append(testFilename,testFileFormat));
            end
        end
    end
end

%% Load data and PF fitting.  Loop over subjects, sessions, filters
nSubjects = length(subjectNameOptions);
 
%% Set all possible filters.
filterOptions = {'A', 'B', 'C', 'D', 'E'};
nFilters = length(filterOptions);

%% Set all possible spatial frequency options
spatialFrequencyOptionsAll = {'3_cpd', '6_cpd', '9_cpd', '12_cpd' '18_cpd'};
maxNSpatialFrequencies = length(spatialFrequencyOptionsAll);

% Fit holding variables.  Not all subjects are run at all sfs,
% so fill with NaNs so it is easier to tell later what we actually
% did.
thresholdFittedRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
slopeFittedRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
thresholdFitted = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
medianThresholdBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
lowThresholdBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
highThresholdBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
slopeFitted = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
medianSlopeBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
lowSlopeBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
highSlopeBoot = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
if (BOOTSTRAP_RAWFITS)
    medianThresholdBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    lowThresholdBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    highThresholdBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    slopeFittedRawRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    medianSlopeBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    lowSlopeBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
    highSlopeBootRaw = NaN*ones(nSubjects,maxNSpatialFrequencies,nFilters);
end

% Set up big lists of what was run.  Want these at full dimension.
subjectBigList = cell(nSubjects,maxNSpatialFrequencies,nFilters);
spatialFrequencyBigList = cell(nSubjects,maxNSpatialFrequencies,nFilters);
filterBigList = cell(nSubjects,maxNSpatialFrequencies,nFilters);
dateBigList = cell(nSubjects,maxNSpatialFrequencies,nFilters);

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
    nSineFreqCyclesPerDegAll(ss) = nSineFreqCyclesPerDeg;
    
    % Set the size of the subplot. Each figure will contain all filters
    % data, so there will be a total of five subplots in one figure.
    sizeSubplot = [2 3];
    
    % Sessions (aka spatial frequency)
    for dd = 1:nSineFreqCyclesPerDeg
        
        % As spatial frequency info is stored in string form, we extract
        % and convert it to double.
        sineFreqCyclesPerDegStr = sineFreqCyclesPerDeg{dd};
        sineFreqCyclesPerDegTemp = sscanf(sineFreqCyclesPerDegStr,'%d');
        
        % Figures for the psychometric function plots
        rawFig = figure; clf;
        set(gcf,'position',[0,0,1920,1080])
        psychoFig = figure; clf;
        set(gcf,'position',[0,0,1920,1080]);

        % Do initial fit for each filter
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
            
            % Load the contrast range data.
            testFilename = GetMostRecentFileName(testFiledir,...
                sprintf('ContrastRange_%s_%d',subjectName,sineFreqCyclesPerDegTemp));
            theContrastData = load(testFilename);
            
            % Extract the threshold from the initial measurements.
            thresholdInitial = theContrastData.preExpDataStruct.thresholdFoundRawLinear;
            
            % Get the tested contrast range here. If there is not the
            % number of contrast points, we just match the size of the
            % array.
            nContrastPoints = 8;
            if ~(length(theContrastData.estDomainValidation)==nContrastPoints)
                theContrastData.estDomainValidation(nContrastPoints) = -10;
            end                
            contrastRangePerSubject(dd,:,ss) = 10.^theContrastData.estDomainValidation;

            % Quest estimator.  Need this here to get data out
            [threshold, para, dataOut] = theData.estimator.thresholdMLE(...
                'thresholdCriterion', thresholdCriterion, 'returnData', true);
            
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

            % Set the contrast levels in linear unit.
            examinedContrastsLinear = 10.^dataOut.examinedContrasts;
            
            % First fit without bootstrapping or plotting to get slopes for each filter
            figure(rawFig);
             if (SUBPLOT)
                subplot(sizeSubplot(1),sizeSubplot(2),ff); hold on;
             end
            if (~BOOTSTRAP_RAWFITS)
                [paramsFittedRaw(:,ff),thresholdFittedRaw(ss,dd,ff), ...
                    ~,~,~, ...
                    ~,~,~,~, ...
                    legendHandles] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
                    'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
                    'newFigureWindow', false, 'pointSize', pointSize, 'axisLog', axisLog,...
                    'questPara', [],'addLegend',false, ...
                    'beta',slopeValList,'nBootstraps',0,'bootConfInterval',bootConfInterval);
            else
                [paramsFittedRaw(:,ff),thresholdFittedRaw(ss,dd,ff), ...
                    medianThresholdBootRaw(ss,dd,ff),lowThresholdBootRaw(ss,dd,ff),highThresholdBootRaw(ss,dd,ff), ...
                    slopeFittedRaw(ss,dd,ff),medianSlopeBootRaw(ss,dd,ff),lowSlopeBootRaw(ss,dd,ff),highSlopeBootRaw(ss,dd,ff), ...
                    legendHandles] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
                    'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
                    'newFigureWindow', false, 'pointSize', pointSize, 'axisLog', axisLog,...
                    'questPara', [],'addLegend',false, ...
                    'beta',slopeValList,'nBootstraps',nBootstraps,'bootConfInterval',bootConfInterval);
            end
            subtitle(sprintf('Raw Fit %d cpd / Filter = %s',sineFreqCyclesPerDegTemp,whichFilter),'fontsize', 12);
            slopeFittedRaw(ss,dd,ff) = paramsFittedRaw(2,ff);

            % Add the entire test image contrast range and chosen ones for
            % the experiment to the plot.
            testContrastMax = max(theContrastData.preExpDataStruct.rawData.testContrast);
            testContrastMin = min(theContrastData.preExpDataStruct.rawData.testContrast);
            nTestContrasts = 30;
            testContrastsLinear = logspace(log10(testContrastMin),log10(testContrastMax),nTestContrasts);
            testContrastsLog = log10(testContrastsLinear);

            % Indicate available and used contrasts
            plot(testContrastsLog,0,'ko','markersize',7);
            plot(theContrastData.estDomainValidation,0,'ko','markersize',7,'markerfacecolor','k');

            % Fussy legend adding
            if (BOOTSTRAP_RAWFITS)
                if (addQuestFit)
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold','BS-Threshold','BS-ConfInt', 'Quest-fit', ...
                        'FontSize', 12, 'location', 'southwest');
                else
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold','BS-Threshold','BS-ConfInt', ...
                        'FontSize', 12, 'location', 'southwest');
                end
            else
                if (addQuestFit)
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold', ...
                        'FontSize', 12, 'location', 'southwest');
                else
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold', ...
                        'FontSize', 12, 'location', 'southwest');
                end
            end

            % Set the range for the x-axis.
            xlim([-3.3 -1]);

            % Force draw
            drawnow;

            % Clear the pointsize for next plot.
            clear pointSize;
        end

        % Finish up and save the raw fit plot
        testFileNameImagesRefine = strrep(theData.describe.testFileNameImages,'-','/');
        testFileNameImagesRefine = strrep(testFileNameImagesRefine,'_','/');

        testFileNameContrastRefine = strrep(theData.describe.testFileNameContrast,'-','/');
        testFileNameContrastRefine = strrep(testFileNameContrastRefine,'_','/');

        main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
        text(0.7,0.4,sprintf('* Subject %s',subjectName),'fontsize',15,'Parent',main);
        text(0.7,0.35,sprintf('* Image file used: %s',testFileNameImagesRefine),'fontsize',15,'Parent',main);
        text(0.7,0.3,sprintf('* Contrast range used (MOA): %s',testFileNameContrastRefine),'fontsize',15,'Parent',main);
 
        % Get the date of experiment.
        testFileNameContrast = theData.describe.testFileNameContrast;
        numExtract = regexp(testFileNameContrast,'\d+','match');
        yearStr = numExtract{3};
        monthStr = numExtract{4};
        dayStr = numExtract{5};
        dateStr = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);

        % Save the plot here if you want.
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
                    sprintf('CS_%s_%d_cpd_%s_raw',subjectName,sineFreqCyclesPerDegTemp,dateStr));
                testFileFormat = '.tiff';
                saveas(gcf,append(testFilename,testFileFormat));
                fprintf('\t Plot has been saved successfully! \n');
                close(gcf);
            end
        end

        % Now fit again using range of slopes determined from first pass.
        % This redoes a lot of things in the loop above that we really
        % shouldn't redo, but it will take a little time to clean up the
        % code.
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

            % Load the contrast range data.
            testFilename = GetMostRecentFileName(testFiledir,...
                sprintf('ContrastRange_%s_%d',subjectName,sineFreqCyclesPerDegTemp));
            theContrastData = load(testFilename);

            % Extract the threshold from the initial measurements.
            thresholdInitial = theContrastData.preExpDataStruct.thresholdFoundRawLinear;

            % Get the tested contrast range here. If there is not the
            % number of contrast points, we just match the size of the
            % array.
            nContrastPoints = 8;
            if ~(length(theContrastData.estDomainValidation)==nContrastPoints)
                theContrastData.estDomainValidation(nContrastPoints) = -10;
            end
            contrastRangePerSubject(dd,:,ss) = 10.^theContrastData.estDomainValidation;

            % Quest estimator.  Need this here to get data out
            [threshold, para, dataOut] = theData.estimator.thresholdMLE(...
                'thresholdCriterion', thresholdCriterion, 'returnData', true);

            % Pull out the data here.
            nTrials = theData.estimator.nRepeat;
            [stimVec, responseVec, structVec] = combineData(theData.estimator);

            % Set the contrast levels in linear unit.
            examinedContrastsLinear = 10.^dataOut.examinedContrasts;

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

            % Quest estimator
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

            % PF fitting here. Use the initial fit to determine slope range
            figure(psychoFig);
            if (SUBPLOT)
                subplot(sizeSubplot(1),sizeSubplot(2),ff); hold on;
            end

            % Figure out slope range based on raw fits
            rawSlopes = squeeze(slopeFittedRaw(ss,dd,:));
            temp = sort(rawSlopes);
            useSlopes = temp(2:nFilters-1);
            minSlopeRange = min(useSlopes)*(10^-slopeRangeLogUnits);
            maxSlopeRange = max(useSlopes)*(10^slopeRangeLogUnits);

            [paramsFitted(:,ff), thresholdFitted(ss,dd,ff), ...
                medianThresholdBoot(ss,dd,ff),lowThresholdBoot(ss,dd,ff),highThresholdBoot(ss,dd,ff), ...
                slopeFitted(ss,dd,ff),medianSlopeBoot(ss,dd,ff),lowSlopeBoot(ss,dd,ff),highSlopeBoot(ss,dd,ff), ...
                legendHandles] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
                'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
                'newFigureWindow', ~SUBPLOT, 'pointSize', pointSize, 'axisLog', axisLog,...
                'questPara', questPara,'addLegend',false, ...
                'beta',linspace(minSlopeRange,maxSlopeRange,nSlopes),'nBootstraps',nBootstraps,'bootConfInterval',bootConfInterval);
            subtitle(sprintf('Second Pass %d cpd / Filter = %s',sineFreqCyclesPerDegTemp,whichFilter),'fontsize', 12);
            slopeFitted(ss,dd,ff) = paramsFitted(2,ff);

            % Add initial threhold to the plot.
            if (addInitialThresholdEst)
                % Plot it on log space if you want.
                if(axisLog)
                    thresholdInitial = log10(thresholdInitial);
                end

                % Plot it here.
                h_high = plot([thresholdInitial(1) thresholdInitial(1)], [0 1], 'b-', 'linewidth',1);
                h_low = plot([thresholdInitial(2) thresholdInitial(2)], [0 1], 'c--', 'linewidth',1);
                plot([thresholdInitial(3) thresholdInitial(3)], [0 1], 'b-', 'linewidth',1);
                plot([thresholdInitial(4) thresholdInitial(4)], [0 1], 'c--', 'linewidth',1);
            end

            % Add the entire test image contrast range and chosen ones for
            % the experiment to the plot.
            testContrastMax = max(theContrastData.preExpDataStruct.rawData.testContrast);
            testContrastMin = min(theContrastData.preExpDataStruct.rawData.testContrast);
            nTestContrasts = 30;
            testContrastsLinear = logspace(log10(testContrastMin),log10(testContrastMax),nTestContrasts);
            testContrastsLog = log10(testContrastsLinear);

            % Indicate available and used contrasts
            plot(testContrastsLog,0,'ko','markersize',7);
            plot(theContrastData.estDomainValidation,0,'ko','markersize',7,'markerfacecolor','k');

            % Fussy legend adding
            if (addInitialThresholdEst)
                if (addQuestFit)
                    legend([legendHandles h_high h_low], 'Data','PF-fit','PF-Threshold','BS-Threshold','BS-ConfInt', 'Quest-fit', ...
                        'Adj From High','Adj From Low', ...
                        'FontSize', 12, 'location', 'southwest');
                else
                    legend([legendHandles h_high h_low], 'Data','PF-fit','PF-Threshold','BS-Threshold','BS-ConfInt', ...
                        'Adj From High','Adj From Low', ...
                        'FontSize', 12, 'location', 'southwest');
                end
            else
                if (addQuestFit)
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold','BS-Threshold','BS-ConfInt', 'Quest-fit', ...
                        'FontSize', 12, 'location', 'southwest');
                else
                    legend(legendHandles, 'Data','PF-Fit','PF-Threshold','BS-Threshold','BS-ConfInt', ...
                        'FontSize', 12, 'location', 'southwest');
                end
            end

            % Set the range for the x-axis.
            xlim([-3.3 -1]);

            % Force draw
            drawnow;

            % Clear the pointsize for next plot.
            clear pointSize;

            % Save up subject, spatial frequency, and filter
            subjectBigList{ss,dd,ff} = subjectName;
            spatialFrequencyBigList{ss,dd,ff} = sineFreqCyclesPerDegTemp;
            filterBigList{ss,dd,ff} = filterOptions{ff};
            dateBigList{ss,dd,ff} = dateStr;
        end

        % Add some text info in the figure.
        %
        % Refine the text before adding it to the figure.
        testFileNameImagesRefine = strrep(theData.describe.testFileNameImages,'-','/');
        testFileNameImagesRefine = strrep(testFileNameImagesRefine,'_','/');

        testFileNameContrastRefine = strrep(theData.describe.testFileNameContrast,'-','/');
        testFileNameContrastRefine = strrep(testFileNameContrastRefine,'_','/');

        main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
        text(0.7,0.4,sprintf('* Subject %s',subjectName),'fontsize',15,'Parent',main);
        text(0.7,0.35,sprintf('* Image file used: %s',testFileNameImagesRefine),'fontsize',15,'Parent',main);
        text(0.7,0.3,sprintf('* Contrast range used (MOA): %s',testFileNameContrastRefine),'fontsize',15,'Parent',main);

        % Print out the progress.
        fprintf('\t Fitting progress - Subject (%d/%d) / Spatial frequency (%d/%d) \n',ss, nSubjects, dd, nSineFreqCyclesPerDeg);

        % Get the date of experiment.
        testFileNameContrast = theData.describe.testFileNameContrast;
        numExtract = regexp(testFileNameContrast,'\d+','match');
        yearStr = numExtract{3};
        monthStr = numExtract{4};
        dayStr = numExtract{5};
        dateStr = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);

        % Save the plot here if you want.
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
                    sprintf('CS_%s_%d_cpd_%s',subjectName,sineFreqCyclesPerDegTemp,dateStr));
                testFileFormat = '.tiff';
                saveas(gcf,append(testFilename,testFileFormat));
                fprintf('\t Plot has been saved successfully! \n');
                close(gcf);
            end
        end

        % Make a table that contains test image profile.
        %
        % Get the number of spatial frequency tested for the previous
        % subject so that we can make a right gap between two adjacent
        % subjects in the table.
        if (ss >= 2)
            nSineFreqCyclesPerDegOneBefore = nSineFreqCyclesPerDegAll(1:ss-1);
            numSpace = sum(nSineFreqCyclesPerDegOneBefore);
        else
            numSpace = 0;
        end

        % Collect the data for the table.
        Date{dd+numSpace,:} = dateStr;
        Subject{dd+numSpace,:} = subjectName;
        SpatialFrequency(dd+numSpace,:) = sineFreqCyclesPerDegTemp;
        
        % Get target primary contrast. We will read it from the test image file
        % we used.
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');   
            testFilenameImage = theData.describe.testFileNameImages;
            fileFormat = '.mat';
            theImageData = load(fullfile(testFiledir,append(testFilenameImage,fileFormat)));
        else
            error('Cannot find data file');
        end
        PrimaryContrast(dd+numSpace) = theImageData.colorDirectionParams.targetScreenPrimaryContrasts(1); 
        clear theImageData;
        
        TestImageContrastMax(dd+numSpace,:) = max(theContrastData.preExpDataStruct.rawData.testContrast);
        
        if strcmp(monthStr,'10')
            RuleMOA{dd+numSpace,:} = '-0.5to+0.3';
        elseif (strcmp(monthStr,'11') & any(strcmp(dayStr,append('0',string([1:1:6])))))
            RuleMOA{dd+numSpace,:} = '-0.5to+0.3';
        else
            RuleMOA{dd+numSpace,:} = '-0.4to+0.4';
        end

        TestContrastNominalMax(dd+numSpace,:) = round(theContrastData.preExpDataStruct.estDomainValidationNominalLinear(end),4);
        TestContrastNominalMin(dd+numSpace,:) = round(theContrastData.preExpDataStruct.estDomainValidationNominalLinear(1),4);
        TestContrasts{dd+numSpace,:} = round(10.^theContrastData.estDomainValidation,4);
    end
end

% Record test image info to excel file.
if (FITALLATONCE)
    if (RECORDTESTIMAGEPROFILE)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
            testFilename = fullfile(testFiledir,'TestImageProfile.xlsx');
        end
        tableImageProfile = table(Date,Subject,SpatialFrequency,PrimaryContrast,TestImageContrastMax,...
            RuleMOA,TestContrastNominalMin,TestContrastNominalMax,TestContrasts);

        % Write a table to the excel file.
        sheet = 1;
        range = 'B2';
        writetable(tableImageProfile,testFilename,'Sheet',sheet,'Range',range);
    end
    fprintf('\t Test image profile has been successfully recorded! \n');
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

% Find cases where the raw threshold is not within the confidence interval
% for the second pass threshold.
indexBigDiff = find(thresholdFittedRaw < lowThresholdBoot | thresholdFittedRaw > highThresholdBoot);

% Summary plots comparing raw and final fits
index = ~isnan(thresholdFittedRaw);
threshFig = figure; clf;
subplot(1,2,1);  hold on
plot(log10(thresholdFittedRaw(index)),log10(thresholdFitted(index)),'bo','MarkerSize',9,'MarkerFaceColor','b');
h = errorbarY(log10(thresholdFittedRaw(index)),log10(medianThresholdBoot(index)), ...
    log10(medianThresholdBoot(index))-log10(lowThresholdBoot(index)),log10(highThresholdBoot(index))-log10(medianThresholdBoot(index)),'go');
set(h,'Color','g'); set(h,'MarkerSize',8); set(h,'MarkerFaceColor','g');
plot(log10(thresholdFittedRaw(indexBigDiff)),log10(thresholdFitted(indexBigDiff)),'ro','MarkerSize',7,'MarkerFaceColor','r');
plot([-3 12],[-3 12],'k');
xlim([-3 8]); ylim([-3 1]);
xlabel('First pass log10 threshold');
ylabel('Second pass log10 threshold');
subplot(1,2,2); hold on
plot(log10(thresholdFittedRaw(index)),log10(thresholdFitted(index)),'bo','MarkerSize',9,'MarkerFaceColor','b');
h = errorbarY(log10(thresholdFittedRaw(index)),log10(medianThresholdBoot(index)), ...
    log10(medianThresholdBoot(index))-log10(lowThresholdBoot(index)),log10(highThresholdBoot(index))-log10(medianThresholdBoot(index)),'go');
set(h,'Color','g'); set(h,'MarkerSize',8); set(h,'MarkerFaceColor','g');
plot(log10(thresholdFittedRaw(indexBigDiff)),log10(thresholdFitted(indexBigDiff)),'ro','MarkerSize',7,'MarkerFaceColor','r');
plot([-3 12],[-3 12],'k');
xlim([-2.5 -0.5]); ylim([-2.5 -0.5]);
xlabel('First pass log10 threshold');
ylabel('Second pass log10 threshold');
saveas(gcf,fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),append('SecondVsFirstPassLogThresholds',testFileFormat)));

slopeFig = figure; clf; hold on
plot(slopeFittedRaw(~isnan(slopeFittedRaw)),slopeFitted(~isnan(slopeFitted)),'ro','MarkerSize',8,'MarkerFaceColor','r');
plot([0 16],[0 16],'k');
xlim([0 16]); ylim([0 16]);
xlabel('First pass slope');
ylabel('Second pass slope');
axis('square');
saveas(gcf,fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),append('SecondVsFirstPassSlopes',testFileFormat)));

% Confidence interval figure;
confIntervalsLog = log10(highThresholdBoot) - log10(lowThresholdBoot);
confIntervalsLin = highThresholdBoot - lowThresholdBoot;
confHist = figure; clf; hold on;
hist(confIntervalsLog(~isnan(confIntervalsLog)),50);
xlabel('Threshold Confidence Interval Magnitude (Log10)');
saveas(gcf,fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),append('ConfidenceIntervalHist',testFileFormat)));

% Bootstrap on raw?
if (BOOTSTRAP_RAWFITS)
    confIntervalsRawLog = log10(highThresholdBootRaw) - log10(lowThresholdBootRaw);
    confIntervalsRawLin = highThresholdBootRaw - lowThresholdBootRaw;
end

% Let's write out the confidence intervals in a sorted order
clear writeCellArray
[~,sortIndex] = sort(confIntervalsLog(:),'descend');
writeCellArray{1,1} = 'Subject';
writeCellArray{1,2} = 'Spatial Freq';
writeCellArray{1,3} = 'Filter';
writeCellArray{1,4} = 'Date';
writeCellArray{1,5} = 'Threshold';
writeCellArray{1,6} = 'Log10 Threshold';
writeCellArray{1,7} = 'Slope';
writeCellArray{1,8} = sprintf('Conf (%d%%)',round(100*bootConfInterval));
writeCellArray{1,9} = sprintf('Log10 Conf (%d%%)',round(100*bootConfInterval));
writeCellArray{1,10} = 'Threshold Raw';
writeCellArray{1,11} = 'Log10 Threshold Raw';
writeCellArray{1,12} = 'Slope Raw';
if (BOOTSTRAP_RAWFITS)
    writeCellArray{1,13} = sprintf('Conf Raw (%d%%)',round(100*bootConfInterval));
    writeCellArray{1,14} = sprintf('Log10 Conf Raw (%d%%)',round(100*bootConfInterval));
end

cellIndex = 2;
for ii = 1:length(sortIndex)
    if (~isnan(confIntervalsLog(sortIndex(ii))))
        writeCellArray{cellIndex,1} = subjectBigList{sortIndex(ii)};
        writeCellArray{cellIndex,2} = spatialFrequencyBigList{sortIndex(ii)};
        writeCellArray{cellIndex,3} = filterBigList{sortIndex(ii)};
        writeCellArray{cellIndex,4} = dateBigStr{sortIndex(ii)};
        writeCellArray{cellIndex,5} = thresholdFitted(sortIndex(ii));
        writeCellArray{cellIndex,6} = log10(thresholdFitted(sortIndex(ii)));
        writeCellArray{cellIndex,7} = slopeFitted(sortIndex(ii));
        writeCellArray{cellIndex,8} = confIntervalsLin(sortIndex(ii));
        writeCellArray{cellIndex,9} = confIntervalsLog(sortIndex(ii));
        writeCellArray{cellIndex,10} = thresholdFittedRaw(sortIndex(ii));
        writeCellArray{cellIndex,11} = log10(thresholdFittedRaw(sortIndex(ii)));
        writeCellArray{cellIndex,12} = slopeFittedRaw(sortIndex(ii));
        if (BOOTSTRAP_RAWFITS)
            writeCellArray{cellIndex,13} = confIntervalsRawLin(sortIndex(ii));
            writeCellArray{cellIndex,14} = confIntervalsRawLog(sortIndex(ii));
        end
        cellIndex = cellIndex+1;
    end
end
writecell(writeCellArray,fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),append('DataSummary','.xlsx')), ...
    'WriteMode',"replacefile");

% Save out full run info
close all
save(fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),'CSFAnalysisOutput'));

