%% SACC_FitCSF
%
% This is to fit CSF curve for SACC project.
%
% See also:
%    asymmetricParabolicFunc

% History:
%    01/13/23   smo    - Started on it.
%    01/19/23   smo    - First fitting CSF with the function and it seems
%                        pretty good. Will be elaborated.
%    02/03/23   smo    - Added a feature to use bootstrapped values to find
%                        a smoothing parameter for smooth spline curve.

%% Initialize.
clear; close all;

%% Set plotting options.
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
thresholdFittedBootRaw = theData.thresholdFittedBootRaw;
medianThresholdBootRaw = theData.medianThresholdBootRaw;
lowThresholdBootRaw = theData.lowThresholdBootRaw;
highThresholdBootRaw = theData.highThresholdBootRaw;

%% Fit CSF.
%
% Get the number of available subjects and filters.
nSubjects = length(theData.subjectNameOptions);
nFilters = length(theData.filterOptions);

% Fitting happens here one by one per subject.
for ss = 1:nSubjects
    % Set a target subject.
    subjectName = theData.subjectNameOptions{ss};
    
    % Set available spatial frequency data for the subject.
    sineFreqCyclesPerDeg = theData.spatialFrequencyOptions(:,ss);
    
    % Here we remove empty cells. It shows a cell empty when a subject does
    % not have all possible spatial frequency data.
    sineFreqCyclesPerDeg = sineFreqCyclesPerDeg(...
        find(~cellfun(@isempty,sineFreqCyclesPerDeg)));
    
    % The number of available data for the subject. It should be 5 to
    % proceed the fitting.
    nSineFreqCyclesPerDeg = length(sineFreqCyclesPerDeg);
    
    % Run fitting only if there are all spatial frequency data.
    if (nSineFreqCyclesPerDeg == 5)
        % Get spatial frequency data in double.
        for dd = 1:nSineFreqCyclesPerDeg
            sineFreqCyclesPerDegNum(dd) = sscanf(sineFreqCyclesPerDeg{dd},'%d');
        end
        
        %% Make a new plot per each subject.
        if (DRAWONEFIGUREPERSUB)
            % Set figure size.
            figureSize = 550;
            
            % Data figure info.
            dataFig = figure; clf; hold on;
            figurePosition = [100 100 figureSize figureSize];
            set(gcf,'position',figurePosition);
            
            % Cross-validation figure info.
            crossFig = figure;
            figurePositionCross = [100+figureSize 100 figureSize figureSize];
            set(gcf,'position',figurePositionCross);
        end
        
        % Here we read out five values of the thresholds (so, five spatial
        % frequency) to fit CSF curve.
        for ff = 1:nFilters
            % Make a new plot per each filter of the subject.
            if (~DRAWONEFIGUREPERSUB)
                % Set figure size.
                figureSize = 550;
                
                % Data figure info.
                dataFig = figure; clf; hold on;
                figurePosition = [100 100 figureSize figureSize];
                set(gcf,'position',figurePosition);
                
                % Cross-validation figure info.
                crossFig = figure;
                figurePositionCross = [100+figureSize 100 figureSize figureSize];
                set(gcf,'position',figurePositionCross);
            end
            
            % Each variable should have the number of 5 entries.
            thresholds = thresholdFittedRaw(ss,:,ff);
            thresholdsBoot = thresholdFittedBootRaw(ss,:,ff,:);
            medianThresholdsBoot = medianThresholdBootRaw(ss,:,ff);
            lowThresholdBoot = lowThresholdBootRaw(ss,:,ff);
            highThresholdBoot = highThresholdBootRaw(ss,:,ff);
            
            % Convert NaN to 0 here.
            for tt = 1:length(thresholds)
                if isnan(thresholds(tt))
                    thresholds(tt) = 0;
                end
            end
            for tt = 1:length(medianThresholdsBoot)
                if isnan(medianThresholdsBoot(tt))
                    medianThresholdsBoot(tt) = 0;
                end
            end
            
            %% Calculate sensitivity.
            sensitivityRawLinear = 1./thresholds;
            sensitivityRawLog = log10(sensitivityRawLinear);
            
            sensitivityMedianBootLinear = 1./medianThresholdsBoot;
            sensitivityMedianBootLog = log10(sensitivityMedianBootLinear);
            
            % All bootstrapped values.
            sensitivityBootLinear = 1./squeeze(thresholdsBoot);
            sensitivityBootLog = log10(sensitivityBootLinear);
            
            % For calculation of confindence interval from bootstrap,
            % (low) threshold becomes (high) sensitivity, and vice
            % versa.
            sensitivityBootHighLinear = 1./lowThresholdBoot;
            sensitivityBootHighLog = log10(sensitivityBootHighLinear);
            
            sensitivityBootLowLinear = 1./highThresholdBoot;
            sensitivityBootLowLog = log10(sensitivityBootLowLinear);
            
            % Calculate spatial frequency in log space.
            sineFreqCyclesPerDegLog = log10(sineFreqCyclesPerDegNum);
            
            %% Sort each array in an ascending order of spatial frequency.
            [sineFreqCyclesPerDegNumSorted I] = sort(sineFreqCyclesPerDegNum,'ascend');
            sineFreqCyclesPerDegLogSorted = sineFreqCyclesPerDegLog(I);
            
            % Linear sorted according to spatial frequency.
            sensitivityRawLinearSorted = sensitivityRawLinear(I);
            sensitivityMedianBootLinearSorted = sensitivityMedianBootLinear(I);
            sensitivityBootHighLinearSorted = sensitivityBootHighLinear(I);
            sensitivityBootLowLinearSorted = sensitivityBootLowLinear(I);
            sensitivityBootLinearSorted = sensitivityBootLinear(I,:);
            
            % Log sorted according to spatial frequency.
            sensitivityRawLogSorted = sensitivityRawLog(I);
            sensitivityMedianBootLogSorted = sensitivityMedianBootLog(I);
            sensitivityBootHighLogSorted = sensitivityBootHighLog(I);
            sensitivityBootLowLogSorted = sensitivityBootLowLog(I);
            sensitivityBootLogSorted = sensitivityBootLog(I,:);
            
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
            
            %% Fitting method 2) Smooth spline method.
            %
            % Load all bootstrapped values.
            if exist('sensitivityBootLinearSorted')
                myCSValsBoot = sensitivityBootLinearSorted';
                nBootPoints = length(sensitivityBootLinearSorted);
                
            else
                % If there is no bootstrapped data, we make the number of
                % points between high and low values.
                myCSValsBootLow = sensitivityBootLowLinearSorted;
                myCSValsBootHigh = sensitivityBootHighLinearSorted;
                nBootPoints = 100;
                
                for bb = 1:length(sensitivityBootLowLinearSorted)
                    myCSValsBoot(:,bb) = linspace(myCSValsBootLow(bb), myCSValsBootHigh(bb), nBootPoints);
                end
            end
            
            % Set the smoothing paramter. You can use bootstrapped values
            % if you want.
            crossValidate = true;
            crossValBoot = true;
            
            if (crossValidate)
                % Set the number of points for plotting the results.
                nSmoothPoints = 100;
                
                % Set the smoothing param searching options.
                nSmoothingParams = 100;
                minSmoothingParam = 0;
                maxSmoothingParam = 0.3;
                crossSmoothingParams = linspace(minSmoothingParam,maxSmoothingParam,nSmoothingParams);
                
                % Make a loop for testing all set smoothing params.
                for sss = 1:length(crossSmoothingParams)
                    smoothCrossError(sss) = 0;
                    for cc = 1:length(mySFVals)
                        % Get the values for the training set.
                        crossIndex = setdiff(1:length(mySFVals),cc);
                        crossFitSFVals = mySFVals(crossIndex);
                        crossFitCSVals = myCSVals(crossIndex);
                        
                        % Fit curve with the training set.
                        smoothFitCross = fit(crossFitSFVals',crossFitCSVals','smoothingspline','SmoothingParam',crossSmoothingParams(sss));
                        
                        % Get the predicted result of testing value.
                        smoothDataPredsCross = feval(smoothFitCross,mySFVals(cc)');
                        
                        % Calculate the error.
                        if (crossValBoot)
                            for bb = 1:nBootPoints
                                smoothCrossError(sss) = smoothCrossError(sss) + norm(myWs(cc)' .* (myCSValsBoot(bb,cc)' - smoothDataPredsCross));
                            end
                        else
                            smoothCrossError(sss) = smoothCrossError(sss) + norm(myWs(cc)' .* (myCSVals(cc)' - smoothDataPredsCross));
                        end
                    end
                end
                
                % Set the smoothing params that has the smallest error.
                [~,index] = min(smoothCrossError);
                smoothingParam = crossSmoothingParams(index);
                
                % Plot cross validation error.
                figure(crossFig); hold on;
                plot(crossSmoothingParams,smoothCrossError,'ko','MarkerSize',6);
                plot(smoothingParam,smoothCrossError(index),'co','MarkerSize',8,'Markerfacecolor','c','Markeredgecolor','k');
                xlabel('Smoothing parameter','fontsize',15);
                ylabel('Cross-validation errors','fontsize',15);
                title('Cross-validation error accoring to smoothing parameter','fontsize',15);
                legend('All params', 'Set param','fontsize',12);
                xlim([minSmoothingParam maxSmoothingParam]);
                
                % Plot the PF figure.
                figure(dataFig);
            else
                % Set the smoothing param as fixed value otherwise.
                smoothingParam = 0.1;
            end
            
            % Get the values for plotting smooth spline fitting curve.
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
            colorOptionsCSF2 = {'k--','r--','g--','b--','c--'};
            colorOptionsCI = {'k','r','g','b','c'};
            
            % Plot raw data.
            plot(sineFreqCyclesPerDegNumSorted, sensitivityRawLinearSorted, colorOptionsRaw{ff},'markersize',20);
            plot(sineFreqCyclesPerDegNumSorted, sensitivityMedianBootLinearSorted, colorOptionsBoot{ff},'markersize',20);
            
            % Plot the CSF fitting with smoothing spline (Method 2).
            if (DRAWONEFIGUREPERSUB)
                plot(smoothPlotSFVals,smoothPlotPreds,colorOptionsCSF2{ff},'LineWidth',4);
            else
                plot(smoothPlotSFVals,smoothPlotPreds,'c-','LineWidth',4);
            end
            
            % Plot the CSF fitting results (Method 1).
            SF_CSF_start = 3;
            SF_CSF_end = 18;
            nPointsCSF = 100;
            sineFreqCyclesPerDegNumCSF = linspace(SF_CSF_start,SF_CSF_end,nPointsCSF);
            sensitivityCSFLinear = asymmetricParabolicFunc(p_optimized, sineFreqCyclesPerDegNumCSF);
            plot(sineFreqCyclesPerDegNumCSF, sensitivityCSFLinear, colorOptionsCSF{ff}, 'linewidth', 2);
            
            % Plot Confidence Interval acquired from Bootstrap.
            errorNeg = abs(sensitivityMedianBootLinearSorted - sensitivityBootLowLinearSorted);
            errorPos = abs(sensitivityBootHighLinearSorted - sensitivityMedianBootLinearSorted);
            errorbar(sineFreqCyclesPerDegNumSorted, sensitivityMedianBootLinearSorted, ...
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
                    close all;
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
                close all;
            end
        end
    end
end
