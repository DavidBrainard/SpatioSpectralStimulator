%% SACC_FitCSF
%
% This is to fit CSF curve for SACC project.
%
% See also:
%    asymmetricParabolicFunc

% History:
%    1/13/23   smo    - Started on it.
%    1/19/23   smo    - First fitting CSF with the function and it seems
%                       pretty good. Will be elaborated.

%% Initialize.
clear; close all;

%% Set the options.
DRAWONEFIGUREPERSUB = false;
WAITFORKEYTODRAW = true;

%% Load the data.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFilename = fullfile(testFiledir,'CSFAnalysisOutput');
    theData = load(testFilename);
    
    % Close the plots if any pops up.
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
nFilters = length(theData.filterOptions);

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
        
        %% Make a new plot per each subject.
        if (DRAWONEFIGUREPERSUB)
            figure; clf; hold on;
            figureSize = 550;
            figurePosition = [100 100 figureSize figureSize];
            set(gcf,'position',figurePosition);
        end
        
        % Here we read out five values of the thresholds (so, five spatial
        % frequency) to fit CSF curve.
        for ff = 1:nFilters
            
            % Make a new plot per each filter of the subject.
            if (~DRAWONEFIGUREPERSUB)
                dataFig = figure; clf; hold on;
                figureSize = 550;
                figurePosition = [100 100 figureSize figureSize];
                set(gcf,'position',figurePosition);
                crossFig = figure;
                figure(dataFig);
            end
            
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
            
            % Linear sorted according to spatial frequency.
            sensitivityRawLinearSorted = sensitivityRawLinear(I);
            sensitivityBootLinearSorted = sensitivityBootLinear(I);
            sensitivityBootHighLinearSorted = sensitivityBootHighLinear(I);
            sensitivityBootLowLinearSorted = sensitivityBootLowLinear(I);
            
            % Log sorted according to spatial frequency.
            sensitivityRawLogSorted = sensitivityRawLog(I);
            sensitivityBootLogSorted = sensitivityBootLog(I);
            sensitivityBootHighLogSorted = sensitivityBootHighLog(I);
            sensitivityBootLowLogSorted = sensitivityBootLowLog(I);
            
            %% Fitting method 1) Asymmetric parabolic function.
            %
            % Set variables to fit CSF.
            mySFVals = sineFreqCyclesPerDegNumSorted;
            myCSVals = sensitivityRawLinearSorted;
            myWs = 1./(sensitivityBootHighLinearSorted-sensitivityBootLowLinearSorted);
            
            % Set parameters for optimazation of the parameter p.
            p0 = [10 0.5 0.1 0.1];
            p_lowerBound = [10 0.5 0.1 0.1];
            p_higherBound = [400 18 10 10];
            A = []; % Set numbers for the condition of A*x <= b*x0 / x0 is the initial point
            b = [];
            Aeq = []; % Matrix for linear equality constraints
            beq = []; % Vector for linear equality constraints
            options = optimset('fmincon');
            
            % Optimize p here.
            p_optimized = fmincon(@(p_unknown) norm(myWs .* (myCSVals - asymmetricParabolicFunc(p_unknown, mySFVals))), ...
                p0, A, b, Aeq, beq, p_lowerBound, p_higherBound, [], options);
            
            %% Fitting method 2) A smooth spline method.
            %
            % It would be helpful to have access to the bootstrapped values
            % for each SF right here.
            %
            
            % See whether fit can take the error bars into account
            nSmoothPoints = 100;
            crossValidate = false;
            if (crossValidate)
                nSmoothingParams = 100;
                maxSmoothingParam = 0.3;
                crossSmoothingParams = linspace(0,maxSmoothingParam,nSmoothingParams);
                for sss = 1:length(crossSmoothingParams)
                    smoothCrossError(sss) = 0;
                    for cc = 1:length(mySFVals)
                        crossIndex = setdiff(1:length(mySFVals),cc);
                        crossFitSFVals = mySFVals(crossIndex);
                        crossFitCSVals = myCSVals(crossIndex);
                        smoothFitCross = fit(crossFitSFVals',crossFitCSVals','smoothingspline','SmoothingParam',crossSmoothingParams(sss));
                        smoothDataPredsCross = feval(smoothFitCross,myCSVals(cc)');
                        smoothCrossError(sss) = smoothCrossError(sss) + norm(myWs(cc)' .* (myCSVals(cc)' - smoothDataPredsCross));
                    end
                end
                figure(crossFig); plot(crossSmoothingParams,smoothCrossError,'ro','MarkerSize',6);
                figure(dataFig);
                [~,index] = min(smoothCrossError);
                smoothingParam = crossSmoothingParams(index)
            else
                smoothingParam = 0.1;
            end
            
            smoothFit = fit(mySFVals',myCSVals','smoothingspline','SmoothingParam',smoothingParam);
            smoothPlotSFVals = linspace(min(mySFVals),max(mySFVals),nSmoothPoints)';
            smoothPlotPreds = feval(smoothFit,smoothPlotSFVals);
            
            %% Plot it here.
            %
            % Different marker/line options for raw fit and bootstrap
            % fit.
            colorOptionsRaw = {'k.','r.','g.','b.','c.'};
            colorOptionsBoot = {'k+','r+','g+','b+','c+'};
            colorOptionsCSF = {'k-','r-','g-','b-','c-'};
            colorOptionsCI = {'k','r','g','b','c'};
            
            % Plot raw data.
            plot(sineFreqCyclesPerDegNumSorted, sensitivityRawLinearSorted, colorOptionsRaw{ff},'markersize',20);
            plot(sineFreqCyclesPerDegNumSorted, sensitivityBootLinearSorted, colorOptionsBoot{ff},'markersize',20);
            
            plot(smoothPlotSFVals,smoothPlotPreds,'c-','LineWidth',4);
            
            % Plot the CSF fitting results.
            SF_CSF_start = 3;
            SF_CSF_end = 18;
            nPointsCSF = 100;
            sineFreqCyclesPerDegNumCSF = linspace(SF_CSF_start,SF_CSF_end,nPointsCSF);
            sensitivityCSFLinear = asymmetricParabolicFunc(p_optimized, sineFreqCyclesPerDegNumCSF);
            plot(sineFreqCyclesPerDegNumCSF, sensitivityCSFLinear, colorOptionsCSF{ff}, 'linewidth', 2);
            
            % Plot Confidence Interval acquired from Bootstrap.
            errorNeg = abs(sensitivityBootLinearSorted - sensitivityBootLowLinearSorted);
            errorPos = abs(sensitivityBootHighLinearSorted - sensitivityBootLinearSorted);
            errorbar(sineFreqCyclesPerDegNumSorted, sensitivityBootLinearSorted, ...
                errorNeg, errorPos, colorOptionsCI{ff});
            
            % Add details per each plot of the subject.
            if (~DRAWONEFIGUREPERSUB)
                xlabel('Spatial Frequency (cpd)','fontsize',15);
                ylabel('Contrast Sensitivity (Linear unit)','fontsize',15);
                
                xticks(sineFreqCyclesPerDegNumSorted);
                xticklabels(sineFreqCyclesPerDegNumSorted);
                
                yaxisRangeLinear = [0:50:300];
                ylim([0 300]);
                yticks(yaxisRangeLinear);
                yticklabels(yaxisRangeLinear);
                
                title(sprintf('CSF curve - Sub %s',subjectName),'fontsize',15);
                
                % Add legend.
                f = flip(get(gca, 'Children'));
                
                numSpaceLegend = 4;
                idxLegendRaw = linspace(1, 1+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
                idxLegendBoot = linspace(2, 2+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
                
                filterWls = {'392 nm' '417 nm' '437 nm' '456 nm' '476 nm'};
                for ll = 1:nSineFreqCyclesPerDeg
                    contentLegendRaw{ll} = sprintf('%s (%s) -PF',theData.filterOptions{ll}, filterWls{ll});
                    contentLegendBoot{ll} = sprintf('%s (%s) -Boot',theData.filterOptions{ll}, filterWls{ll});
                end
                
                % Add legend when drawing one figure per each filter.
                legend(append(filterWls{ff},'-PF'),append(filterWls{ff},'-Boot'),...
                    'fontsize', 13, 'location', 'northeast');
            end
            
            % Key stroke to start measurement.
            if (~DRAWONEFIGUREPERSUB)
                if (WAITFORKEYTODRAW)
                    disp('Press a key to draw next plot!');
                    pause;
                    close;
                end
            end
        end
        
        %% Add details per each plot of the subject.
        if (DRAWONEFIGUREPERSUB)
            xlabel('Spatial Frequency (cpd)','fontsize',15);
            ylabel('Contrast Sensitivity (Linear unit)','fontsize',15);
            
            xticks(sineFreqCyclesPerDegNumSorted);
            xticklabels(sineFreqCyclesPerDegNumSorted);
            
            yaxisRangeLinear = [0:50:300];
            ylim([0 300]);
            yticks(yaxisRangeLinear);
            yticklabels(yaxisRangeLinear);
            
            title(sprintf('CSF curve - Sub %s',subjectName),'fontsize',15);
            
            % Add legend.
            f = flip(get(gca, 'Children'));
            
            numSpaceLegend = 4;
            idxLegendRaw = linspace(1, 1+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            idxLegendBoot = linspace(2, 2+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            
            filterWls = {'392 nm' '417 nm' '437 nm' '456 nm' '476 nm'};
            for ll = 1:nSineFreqCyclesPerDeg
                contentLegendRaw{ll} = sprintf('%s (%s) -PF',theData.filterOptions{ll}, filterWls{ll});
                contentLegendBoot{ll} = sprintf('%s (%s) -Boot',theData.filterOptions{ll}, filterWls{ll});
            end
            
            % Add legend when drawing one figure per each subject.
            legend(f([idxLegendRaw idxLegendBoot]), [contentLegendRaw contentLegendBoot], ...
                'fontsize', 13, 'location', 'northeast');
            
            % Key stroke to start measurement.
            if (WAITFORKEYTODRAW)
                disp('Press a key to draw next plot!');
                pause;
                close;
            end
        end
    end
end
