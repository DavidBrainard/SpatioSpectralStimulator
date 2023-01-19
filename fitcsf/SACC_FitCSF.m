%% SACC_FitCSF
%
% This is to fit CSF curve for SACC project.
%
% See also:
%    asymmetricParabolicFunc

% History:
%    1/13/23   smo    - Started on it.

%% Initialize.
clear; close all;

%% Load the data.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFilename = fullfile(testFiledir,'CSFAnalysisOutput');
    theData = load(testFilename);
    
    % Close the plots if popping up any.
    close all;
else
    error('Cannot find the data file!');
end

%% Read out the data.¸
%
% Subject info.
subjectName = theData.subjectNameOptions;

% Get threshold data. Each variable is aligned in [subject, SF, filter].
thresholdFittedRaw = theData.thresholdFittedRaw;
medianThresholdBootRaw = theData.medianThresholdBootRaw;
lowThresholdBootRaw = theData.lowThresholdBootRaw;
highThresholdBootRaw = theData.highThresholdBootRaw;

%% Fit CSF.
nSubjects = length(theData.subjectNameOptions);
nFilters = 5;

for ss = 1:nSubjects
    %% Set target subject and filter. We will fit CSF one by one.
    subjectName = theData.subjectNameOptions{ss};
    
    % Set available spatial frequency data for the subject.
    sineFreqCyclesPerDeg = theData.spatialFrequencyOptions(:,ss);
    
    % Here we remove empty cells. It shows a cell empty when a subject does
    % not have all possible spatial frequency data.
    sineFreqCyclesPerDeg = sineFreqCyclesPerDeg(...
        find(~cellfun(@isempty,sineFreqCyclesPerDeg)));
    
    nSineFreqCyclesPerDeg = length(sineFreqCyclesPerDeg);
    
    % Run fitting only if there are all spatial frequency data.
    if (nSineFreqCyclesPerDeg == 5)
        % Get spatial frequency data in double.
        for dd = 1:nSineFreqCyclesPerDeg
            sineFreqCyclesPerDegNum(dd) = sscanf(sineFreqCyclesPerDeg{dd},'%d');
        end
        
        % Here we read out five values of the thresholds (so, five spatial
        % frequency) to fit CSF curve.
        for ff = 1:nFilters
            % Each variable should have the number of 5 entries.
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
            
            %% Calculate sensitivity.
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
            
            sensitivityRawLinearSorted = sensitivityRawLinear(I);
            sensitivityRawLogSorted = sensitivityRawLog(I);
            sensitivityBootLogSorted = sensitivityBootLog(I);
            sensitivityBootHighLogSorted = sensitivityBootHighLog(I);
            sensitivityBootLowLogSorted = sensitivityBootLowLog(I);
            
            %% Fitting happens here.
            %
            % Set variables to fit CSF.
            mySFVals = sineFreqCyclesPerDegNumSorted;
            myCSVals = sensitivityRawLinearSorted;
            myWs = 1;
            
            % Set parameters for optimazation of the parameter p.
            p0 = [100 3 1 1];
            p_lowerBound = [10 0.5 0.1 0.1];
            p_higherBound = [500 18 10 10];
            A = []; % Set numbers for the condition of A*x <= b*x0 / x0 is the initial point
            b = []; 
            Aeq = []; % Matrix for linear equality constraints
            beq = []; % Vector for linear equality constraints
            options = optimset('fmincon'); 

            % Optimize p here.
            p_optimized = fmincon(@(p_unknown) norm(myWs .* (myCSVals - asymmetricParabolicFunc(p_unknown, mySFVals))), ...
                p0, A, b, Aeq, beq, p_lowerBound, p_higherBound, [], options);
                      
            %% Plot it here.
            %
            % Make a figure per each subject.
            figure; clf; hold on;
            figureSize = 900;
            figurePosition = [1000 500 figureSize figureSize];
            set(gcf,'position',figurePosition);
            
            % Different marker/line options for raw fit and bootstrap
            % fit.
            colorOptionsRaw = {'k.','r.','g.','b.','c.'};
            colorOptionsBoot = {'k+','r+','g+','b+','c+'};
            colorOptionsCI = {'k','r','g','b','c'};
            
            % Plot raw data.
            plot(sineFreqCyclesPerDegLogSorted, sensitivityRawLogSorted, colorOptionsRaw{ff},'markersize',20,'linewidth',2);
            plot(sineFreqCyclesPerDegLogSorted, sensitivityBootLogSorted, colorOptionsBoot{ff},'markersize',20,'linewidth',2);
            
            % Plot the CSF fitting results.
            sineFreqCyclesPerDegNumCSF = linspace(1,20,20);
            sensitivityCSFLinear = asymmetricParabolicFunc(p_optimized, sineFreqCyclesPerDegNumCSF);
            plot(log10(sineFreqCyclesPerDegNumCSF), log10(sensitivityCSFLinear), 'r-');
            
            % Plot Confidence Interval acquired from Bootstrap.
%             errorNeg = abs(sensitivityBootLogSorted - sensitivityBootLowLogSorted);
%             errorPos = abs(sensitivityBootHighLogSorted - sensitivityBootLogSorted);
%             errorbar(sineFreqCyclesPerDegLogSorted, sensitivityBootLogSorted, ...
%                 errorNeg, errorPos, colorOptionsCI{ff});
        end
        
        %% Add details per each plot.
        % Set axis and stuffs.
        xlabel('Spatial Frequency (cpd)','fontsize',15);
        ylabel('Contrast Sensitivity','fontsize',15);
        
        xticks(sineFreqCyclesPerDegLogSorted);
        xticklabels(sineFreqCyclesPerDegNumSorted);
        
        yaxisRangeLinear = [0:20:300];
        ylim(log10([0 300]));
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
            contentLegendRaw{ll} = sprintf('%s (%s) -PF',theData.filterOptions{ll}, filterWls{ll});
            contentLegendBoot{ll} = sprintf('%s (%s) -Boot',theData.filterOptions{ll}, filterWls{ll});
        end
%         legend(f([idxLegendRaw idxLegendBoot]), [contentLegendRaw contentLegendBoot], ...
%             'fontsize', 13, 'location', 'northeast');
    end
end

%% Save the results.
