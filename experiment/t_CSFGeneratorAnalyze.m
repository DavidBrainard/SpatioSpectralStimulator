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
%    03/09/23  smo          - Removed the part checking Quest+ adaptive
%                             mode. It has been saved in the copy of this
%                             code, t_CSFGeneratorAnalyze_OLD.m. Also,
%                             added the part that creates a text file
%                             containing all CS values per each subject

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
FITALLATONCE = false;
SAVETHEPLOT = false;
RECORDTESTIMAGEPROFILE = true;
RECORDTEXTSUMMARYPERSUB = true;

PLOTCSFCURVE = true;
SAVECSFCURVE = true;

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
                if ~strcmp(dataList(idxData).name,'FlickerPhotom')
                    spatialFrequencyOptions{dd,ss} = dataList(idxData).name;
                end
            end
        end
    end
else
    % Load single subject to fit one by one.
    subjectNameOptions = {'002'};
    spatialFrequencyOptions = {'3'; '6'; '9'; '12'; '18'};
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
        
        % Collect spatial frequency in double format. This will be used
        % when plotting CSF curve when all data is available.
        sineFreqCyclesPerDegNum(dd) = sineFreqCyclesPerDegTemp;
        
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
                theContrastData.estDomainValidation(nContrastPoints) = NaN;
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
                    ~,~,~,~, ...
                    ~,~,~,~,~,~, ...
                    legendHandles] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
                    'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE,'paramsFree', paramsFree, ...
                    'newFigureWindow', false, 'pointSize', pointSize, 'axisLog', axisLog,...
                    'questPara', [],'addLegend',false, ...
                    'beta',slopeValList,'nBootstraps',0,'bootConfInterval',bootConfInterval);
            else
                [paramsFittedRaw(:,ff),thresholdFittedRaw(ss,dd,ff), thresholdFittedBootRaw(ss,dd,ff,:),...
                    medianThresholdBootRaw(ss,dd,ff),lowThresholdBootRaw(ss,dd,ff),highThresholdBootRaw(ss,dd,ff), ...
                    slopeFittedRaw(ss,dd,ff),medianSlopeBootRaw(ss,dd,ff),lowSlopeBootRaw(ss,dd,ff),highSlopeBootRaw(ss,dd,ff), ...
                    thresholdFittedBootCross1Raw(ss,dd,ff,:),thresholdFittedBootCross2Raw(ss,dd,ff,:), ...
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
        
        %% Now fit again using range of slopes determined from first pass.
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
            
            [paramsFitted(:,ff), thresholdFitted(ss,dd,ff), thresholdFittedBoot(ss,dd,ff,:),...
                medianThresholdBoot(ss,dd,ff),lowThresholdBoot(ss,dd,ff),highThresholdBoot(ss,dd,ff), ...
                slopeFitted(ss,dd,ff),medianSlopeBoot(ss,dd,ff),lowSlopeBoot(ss,dd,ff),highSlopeBoot(ss,dd,ff), ...
                thresholdFittedBootCross1(ss,dd,ff,:),thresholdFittedBootCross2(ss,dd,ff,:),...
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
        PrimaryContrast(dd+numSpace,:) = theImageData.colorDirectionParams.targetScreenPrimaryContrasts(1);
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
        
        % Check if the threshold lies within the range of test contrasts.
        % If the threshold was estimated within the test contrast range,
        % print out 'Good', else 'Bad'.
        minTestContrast = min(contrastRangePerSubject(dd,:));
        maxTestContrast = max(contrastRangePerSubject(dd,:));
        if or(any(thresholdFittedRaw(ss,dd,:) < minTestContrast),any(thresholdFittedRaw(ss,dd,:) > maxTestContrast))
            ThresholdEstimate{dd+numSpace,:} = 'Bad';
        else
            ThresholdEstimate{dd+numSpace,:} = 'Good';
        end
    end
    
    %% Plot CSF curve here.
    % This part will run only if all spatial frequency data is avaialable.
    if (PLOTCSFCURVE)
        if (nSineFreqCyclesPerDeg == maxNSpatialFrequencies)
            % Make a new figure to plot CSF.
            figure; clf; hold on;
            figureSize = 550;
            figurePosition = [1000 500 figureSize figureSize];
            set(gcf,'position',figurePosition);
            
            % Export the threshold data.
            for ff = 1:nFilters
                thresholdsRaw = thresholdFittedRaw(ss,:,ff);
                thresholdsBoot = medianThresholdBootRaw(ss,:,ff);
                thresholdBootLow = lowThresholdBootRaw(ss,:,ff);
                thresholdBootHigh = highThresholdBootRaw(ss,:,ff);
                
                % Convert NaN to 0 here.
                for tt = 1:length(thresholdsRaw)
                    if isnan(thresholdsRaw(tt))
                        thresholdsRaw(tt) = 0;
                    end
                end
                for tt = 1:length(thresholdsBoot)
                    if isnan(thresholdsBoot(tt))
                        thresholdsBoot(tt) = 0;
                    end
                end
                
                % Calculate sensitivity.
                sensitivityRawLinear = 1./thresholdsRaw;
                sensitivityRawLog = log10(sensitivityRawLinear);
                
                sensitivityBootLinear = 1./thresholdsBoot;
                sensitivityBootLog = log10(sensitivityBootLinear);
                
                % For calculation of confindence interval from bootstrap,
                % (low) threshold becomes (high) sensitivity, and vice
                % versa.
                sensitivityBootHighLinear = 1./thresholdBootLow;
                sensitivityBootHighLog = log10(sensitivityBootHighLinear);
                
                sensitivityBootLowLinear = 1./thresholdBootHigh;
                sensitivityBootLowLog = log10(sensitivityBootLowLinear);
                
                % Calculate spatial frequency in log space.
                sineFreqCyclesPerDegLog = log10(sineFreqCyclesPerDegNum);
                
                % Sort each array in an ascending order of spatial
                % frequency.
                [sineFreqCyclesPerDegNumSorted I] = sort(sineFreqCyclesPerDegNum,'ascend');
                sineFreqCyclesPerDegLogSorted = sineFreqCyclesPerDegLog(I);
                
                sensitivityRawLogSorted = sensitivityRawLog(I);
                sensitivityBootLogSorted = sensitivityBootLog(I);
                sensitivityBootHighLogSorted = sensitivityBootHighLog(I);
                sensitivityBootLowLogSorted = sensitivityBootLowLog(I);
                
                % Plot it here.
                %
                % Different marker/line options for raw fit and bootstrap
                % fit.
                colorOptionsRaw = {'k.-','r.-','g.-','b.-','c.-'};
                colorOptionsBoot = {'k.--','r.--','g.--','b.--','c.--'};
                colorOptionsCI = {'k','r','g','b','c'};
                
                % Plot CSF curves.
                plot(sineFreqCyclesPerDegLogSorted, sensitivityRawLogSorted, colorOptionsRaw{ff},'markersize',20,'linewidth',2);
                plot(sineFreqCyclesPerDegLogSorted, sensitivityBootLogSorted, colorOptionsBoot{ff},'markersize',20,'linewidth',2);
                
                % Plot Confidence Interval acquired from Bootstrap.
                errorNeg = abs(sensitivityBootLogSorted - sensitivityBootLowLogSorted);
                errorPos = abs(sensitivityBootHighLogSorted - sensitivityBootLogSorted);
                errorbar(sineFreqCyclesPerDegLogSorted, sensitivityBootLogSorted, ...
                    errorNeg, errorPos, colorOptionsCI{ff});
            end
            
            % Set axis and stuffs.
            xlabel('Spatial Frequency (cpd)','fontsize',15);
            ylabel('Contrast Sensitivity','fontsize',15);
            
            xticks(sineFreqCyclesPerDegLogSorted);
            xticklabels(sineFreqCyclesPerDegNumSorted);
            
            yaxisRangeLinear = [0:50:300];
            ylim(log10([1 300]));
            yticks(log10(yaxisRangeLinear));
            yticklabels(yaxisRangeLinear);
            
            title(sprintf('CSF curve - Sub %s',subjectName),'fontsize',15);
            subtitle('Data points are in log unit, labels are in linear unit', 'fontsize',13);
            
            % Add legend.
            f = flip(get(gca, 'Children'));
            
            numSpaceLegend = 3;
            idxLegendRaw = linspace(1, 1+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            idxLegendBoot = linspace(2, 2+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            
            filterWls = {'392 nm' '417 nm' '437 nm' '456 nm' '476 nm'};
            for ll = 1:nSineFreqCyclesPerDeg
                contentLegendRaw{ll} = sprintf('%s (%s) -PF',filterOptions{ll}, filterWls{ll});
                contentLegendBoot{ll} = sprintf('%s (%s) -Boot',filterOptions{ll}, filterWls{ll});
            end
            legend(f([idxLegendRaw idxLegendBoot]), [contentLegendRaw contentLegendBoot], ...
                'fontsize', 13, 'location', 'northeast');
            
            % Save the plot if you want.
            if (SAVECSFCURVE)
                if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),...
                        subjectName,'CSF');
                    
                    % Make folder with subject name if it does not exist.
                    if ~exist(testFiledir, 'dir')
                        mkdir(testFiledir);
                    end
                    
                    % Save the plot.
                    testFilename = fullfile(testFiledir,sprintf('CSF_%s',subjectName));
                    testFileFormat = '.tiff';
                    saveas(gcf,append(testFilename,testFileFormat));
                    fprintf('\t Plot has been saved successfully! \n');
                    close(gcf);
                end
            end
        end
    end
    
    %% Save out summary text file per subject.
    if (RECORDTEXTSUMMARYPERSUB)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),subjectName,'CSF');
            testFilename = fullfile(testFiledir,sprintf('CS_Summary_%s.xlsx',subjectName));
        end
        
        % Sort each data in a single column.
        nCSPerSub = nFilters * nSineFreqCyclesPerDeg;
        
        NumCount_CSSummary = linspace(1,nCSPerSub,nCSPerSub)';
        Subject_CSSummary = reshape(squeeze(subjectBigList(ss,:,:)),nCSPerSub,1);
        Filter_CSSummary = reshape(squeeze(filterBigList(ss,:,:)),nCSPerSub,1);
        SpatialFrequency_CSSummary = reshape(squeeze(spatialFrequencyBigList(ss,:,:)),nCSPerSub,1);
        ThresholdPF_CSSummary = reshape(squeeze(thresholdFittedRaw(ss,:,:)),nCSPerSub,1);
        MedianThresholdBoot_CSSummary = reshape(squeeze(medianThresholdBootRaw(ss,:,:)),nCSPerSub,1);
        BootCILow_CSSummary = reshape(squeeze(lowThresholdBootRaw(ss,:,:)),nCSPerSub,1);
        BootCIHigh_CSSummary = reshape(squeeze(highThresholdBootRaw(ss,:,:)),nCSPerSub,1);
        
        LogSensitivityPF_CSSummary = log10(1./ThresholdPF_CSSummary);
        LogSensitivityMedianBoot_CSSummary = log10(1./MedianThresholdBoot_CSSummary);
        LogSensitivitiyBootCILow_CSSummary = log10(1./BootCILow_CSSummary);
        LogSensitivitiyBootCIHigh_CSSummary = log10(1./BootCIHigh_CSSummary);
        
        % Make a table.
        tableCSSummary = table(NumCount_CSSummary, Subject_CSSummary,Filter_CSSummary, SpatialFrequency_CSSummary, ThresholdPF_CSSummary, ...
            MedianThresholdBoot_CSSummary, BootCILow_CSSummary, BootCIHigh_CSSummary,...
            LogSensitivityPF_CSSummary,LogSensitivityMedianBoot_CSSummary,LogSensitivitiyBootCILow_CSSummary,LogSensitivitiyBootCIHigh_CSSummary);
        
        % Change the variable name as desired.
        tableCSSummary.Properties.VariableNames = {'No', 'Subject', 'Filter', 'SpatialFrequency',...
            'ThresholdPF', 'MedianThresholdBoot', 'BootCILow', 'BootCIHigh',...
            'LogSensitivityPF', 'LogSensitivityMedianThresholdBoot', 'LogSensitivityBootCILow', 'LogSensitivityBootCIHigh'};
        
        % Write a table to the excel file.
        sheet = 1;
        range = 'B2';
        writetable(tableCSSummary,testFilename,'Sheet',sheet,'Range',range);
    end
end

%% Record test image info to excel file.
if (FITALLATONCE)
    if (RECORDTESTIMAGEPROFILE)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'));
            testFilename = fullfile(testFiledir,'TestImageProfile.xlsx');
        end
        tableImageProfile = table(Date,Subject,SpatialFrequency,PrimaryContrast,TestImageContrastMax,...
            RuleMOA,TestContrastNominalMin,TestContrastNominalMax,TestContrasts,ThresholdEstimate);
        
        % Write a table to the excel file.
        sheet = 1;
        range = 'B2';
        writetable(tableImageProfile,testFilename,'Sheet',sheet,'Range',range);
    end
    fprintf('\t Test image profile has been successfully recorded! \n');
end

%% Run this part only when fitting all at once.
%
% This part saves some summary text file and plots, we will not run this
% part when fitting one subject to prevent overwriting the data and figure.
if (FITALLATONCE)
    % Find cases where the raw threshold is not within the confidence
    % interval for the second pass threshold.
    indexBigDiff = find(thresholdFittedRaw < lowThresholdBoot | thresholdFittedRaw > highThresholdBoot);
    
    % Summary plots comparing raw and final fits.
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
            writeCellArray{cellIndex,4} = dateBigList{sortIndex(ii)};
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
    
    % Save out all experiment results.
    writecell(writeCellArray,fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),append('DataSummary','.xlsx')), ...
        'WriteMode',"replacefile");
    
    % Save out full run info.
    save(fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),'CSFAnalysisOutput'));
end

%% Close all.
close all;
