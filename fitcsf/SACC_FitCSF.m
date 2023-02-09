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
%    02/08/23   smo    - Now we can cross-validate with two separate
%                        bootstrapped data when we fit CSF with smooth
%                        spline function.

%% Initialize.
clear; close all;

%% Decide the options to fit CSF and how to plot.
%
% Plotting options.
OneFigurePerSub = false;
WaitForKeyToPlot = true;

% Fitting options.
FitAsymmetricParabolic = false;
FitSmoothSpline = true;

% When fitting Smooth spline, You can choose option among {'crossVal',
% 'crossValBootWithin', 'crossValBootAcross', 'type'}.
if (FitSmoothSpline)
    OptionSearchSmoothParam = 'crossValBootAcross';
end

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
thresholdFittedBootCross1Raw = theData.thresholdFittedBootCross1Raw;
thresholdFittedBootCross2Raw = theData.thresholdFittedBootCross2Raw;

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
        
        %% Make a new plot per each subject if you want.
        if (OneFigurePerSub)
            % Set figure size.
            figureSize = 550;
            
            % Data figure info.
            dataFig = figure; clf; hold on;
            figurePosition = [100 100 figureSize figureSize];
            set(gcf,'position',figurePosition);
            
            % Cross-validation figure info.
            if (FitSmoothSpline)
                crossFig = figure; hold on;
                figurePositionCross = [100+figureSize 100 figureSize figureSize];
                set(gcf,'position',figurePositionCross);
            end
        end
        
        % Here we read out five values of the thresholds (so, five spatial
        % frequency) to fit CSF curve.
        for ff = 1:nFilters
            % Make a new plot per each filter of the subject.
            if (~OneFigurePerSub)
                % Set figure size.
                figureSize = 550;
                
                % Data figure info.
                dataFig = figure; clf; hold on;
                figurePosition = [100 100 figureSize figureSize];
                set(gcf,'position',figurePosition);
                
                % Cross-validation figure info.
                if (FitSmoothSpline)
                    crossFig = figure; hold on;
                    figurePositionCross = [100+figureSize 100 figureSize figureSize];
                    set(gcf,'position',figurePositionCross);
                end
            end
            
            % Read out the variables per each filter.
            thresholds = thresholdFittedRaw(ss,:,ff);
            thresholdsBoot = thresholdFittedBootRaw(ss,:,ff,:);
            medianThresholdsBoot = medianThresholdBootRaw(ss,:,ff);
            lowThresholdBoot = lowThresholdBootRaw(ss,:,ff);
            highThresholdBoot = highThresholdBootRaw(ss,:,ff);
            thresholdsBootCross1 = thresholdFittedBootCross1Raw(ss,:,ff,:);
            thresholdsBootCross2 = thresholdFittedBootCross2Raw(ss,:,ff,:);
            
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
            
            % Additional bootstrapped values for cross-validation.
            sensitivityBootCross1Linear = 1./squeeze(thresholdsBootCross1);
            sensitivityBootCross2Linear = 1./squeeze(thresholdsBootCross2);
            
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
            sensitivityBootCross1LinearSorted = sensitivityBootCross1Linear(I,:);
            sensitivityBootCross2LinearSorted = sensitivityBootCross2Linear(I,:);
            
            % Log sorted according to spatial frequency.
            sensitivityRawLogSorted = sensitivityRawLog(I);
            sensitivityMedianBootLogSorted = sensitivityMedianBootLog(I);
            sensitivityBootHighLogSorted = sensitivityBootHighLog(I);
            sensitivityBootLowLogSorted = sensitivityBootLowLog(I);
            sensitivityBootLogSorted = sensitivityBootLog(I,:);
            
            %% Set variables to fit CSF.
            mySFVals = sineFreqCyclesPerDegNumSorted;
            myCSVals = sensitivityRawLinearSorted;
            myWs = 1./(sensitivityBootHighLinearSorted-sensitivityBootLowLinearSorted);
            
            %% Fitting method 1) Asymmetric parabolic function.
            if (FitAsymmetricParabolic)
                % Set parameters for optimazation of the parameter p.
                p0 = [10 0.5 0.1 0.1];
                p_lowerBound = [10 0.5 0.1 0.1];
                p_higherBound = [400 18 10 10];
                A = [];   % Set numbers for the condition of A*x <= b*x0 / x0 is the initial point
                b = [];
                Aeq = []; % Matrix for linear equality constraints
                beq = []; % Vector for linear equality constraints
                options = optimset('fmincon');
                
                % Optimize p here.
                p_optimized = fmincon(@(p_unknown) norm(myWs .* (myCSVals - asymmetricParabolicFunc(p_unknown, mySFVals))), ...
                    p0, A, b, Aeq, beq, p_lowerBound, p_higherBound, [], options);
            end
            
            %% Fitting method 2) Smooth spline method.
            if (FitSmoothSpline)
                % Load all bootstrapped values.
                myCSValsBoot = sensitivityBootLinearSorted';
                myCSValsCross1 = sensitivityBootCross1LinearSorted';
                myCSValsCross2 = sensitivityBootCross2LinearSorted';
                nBootPoints = size(myCSValsBoot,1);
                
                %% Search smoothing parameter here.
                %
                % Set the number of points for plotting the results.
                nSmoothPoints = 100;
                
                % Set the smoothing param searching options.
                nSmoothingParams = 100;
                minSmoothingParam = 0;
                maxSmoothingParam = 1;
                crossSmoothingParams = linspace(minSmoothingParam,maxSmoothingParam,nSmoothingParams);
                
                % We search smoothing parameter based on the option selected
                % from the above.
                switch OptionSearchSmoothParam
                    case 'crossValBootWithin'
                        % Set the number of bootstrapped points to use.
                        nCrossValBootWithin = 20;
                        
                        for sss = 1:length(crossSmoothingParams)
                            smoothCrossError(sss) = 0;
                            
                            for cc = 1:nCrossValBootWithin
                                % Make the training and validating dataset.
                                % Here we only pick randomly once at the
                                % beginning and keep using the same dataset.
                                if (sss == 1)
                                    for zz = 1:length(mySFVals)
                                        bootCSFDataFit{cc}(zz) = myCSValsBoot(randi(nBootPoints,1,1),zz);
                                        bootCSFDataCross{cc}(zz) = myCSValsBoot(randi(nBootPoints,1,1),zz);
                                    end
                                end
                                
                                % Fit curve with the training set.
                                smoothFitCross = fit(mySFVals',bootCSFDataFit{cc}','smoothingspline','SmoothingParam',crossSmoothingParams(sss));
                                
                                % Get the predicted result of testing value.
                                smoothDataPredsCross = feval(smoothFitCross,mySFVals');
                                
                                % Calculate the error.
                                smoothCrossError(sss) = smoothCrossError(sss) + sum((bootCSFDataCross{cc}' - smoothDataPredsCross).^2);
                            end
                            % Print out the progress.
                            fprintf('Smoothing param testing progress - (%d/%d) \n', sss, length(crossSmoothingParams));
                        end
                        
                    case 'crossValBootAcross'
                        nCrossValBootAcross = 20;
                        for sss = 1:length(crossSmoothingParams)
                            smoothCrossError(sss) = 0;
                            
                            for cc = 1:nCrossValBootAcross
                                if (sss == 1)
                                    for zz = 1:length(mySFVals)
                                        crossIndex = randi(nBootPoints,1,1);
                                        bootCSFDataFit{cc}(zz) = myCSValsCross1(crossIndex,zz);
                                        bootCSFDataCross{cc}(zz) = myCSValsCross2(crossIndex,zz);
                                    end
                                end
                                
                                % Fit curve with the training set.
                                smoothFitCross = fit(mySFVals',bootCSFDataFit{cc}','smoothingspline','SmoothingParam',crossSmoothingParams(sss));
                                
                                % Get the predicted result of testing value.
                                smoothDataPredsCross = feval(smoothFitCross,mySFVals');
                                
                                % Calculate the error.
                                smoothCrossError(sss) = smoothCrossError(sss) + sum((bootCSFDataCross{cc}' - smoothDataPredsCross).^2);
                            end
                            % Print out the progress.
                            fprintf('Smoothing param testing progress - (%d/%d) \n', sss, length(crossSmoothingParams));
                        end
                        
                    case 'crossVal'
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
                                smoothCrossError(sss) = smoothCrossError(sss) + norm(myWs(cc)' .* (myCSVals(cc)' - smoothDataPredsCross));
                            end
                        end
                        
                    case 'type'
                        % Type a number manually.
                        smoothingParam = 0.1;
                end
                
                % Plot the bootstrapped results if we did.
                if any(strcmp(OptionSearchSmoothParam,{'crossValBootWithin','crossValBootAcross','crossVal'}))
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
                end
                
                % Plot the PF figure.
                figure(dataFig);
                
                % Get the values for plotting smooth spline fitting curve.
                smoothFit = fit(mySFVals',myCSVals','smoothingspline','SmoothingParam',smoothingParam);
                smoothPlotSFVals = linspace(min(mySFVals),max(mySFVals),nSmoothPoints)';
                smoothPlotPreds = feval(smoothFit,smoothPlotSFVals);
            end
            
            %% Plot it here.
            %
            % Set marker/line options for the plot.
            colorOptionsRaw = {'k.','r.','g.','b.','c.'};
            colorOptionsBoot = {'k+','r+','g+','b+','c+'};
            colorOptionsCSF = {'k-','r-','g-','b-','c-'};
            colorOptionsCSF2 = {'k--','r--','g--','b--','c--'};
            colorOptionsCI = {'k','r','g','b','c'};
            
            % Raw data.
            plot(sineFreqCyclesPerDegNumSorted, sensitivityRawLinearSorted, colorOptionsRaw{ff},'markersize',20);
            plot(sineFreqCyclesPerDegNumSorted, sensitivityMedianBootLinearSorted, colorOptionsBoot{ff},'markersize',20);
            
            % Confidence Interval.
            errorNeg = abs(sensitivityMedianBootLinearSorted - sensitivityBootLowLinearSorted);
            errorPos = abs(sensitivityBootHighLinearSorted - sensitivityMedianBootLinearSorted);
            e = errorbar(sineFreqCyclesPerDegNumSorted, sensitivityMedianBootLinearSorted, ...
                errorNeg, errorPos, colorOptionsCI{ff});
            e.LineStyle = 'none';
            
            % CSF fitting results with asymmetric parabolic (Method 1).
            if (FitAsymmetricParabolic)
                SF_CSF_start = 3;
                SF_CSF_end = 18;
                nPointsCSF = 100;
                sineFreqCyclesPerDegNumCSF = linspace(SF_CSF_start,SF_CSF_end,nPointsCSF);
                sensitivityCSFLinear = asymmetricParabolicFunc(p_optimized, sineFreqCyclesPerDegNumCSF);
                plot(sineFreqCyclesPerDegNumCSF, sensitivityCSFLinear, colorOptionsCSF{ff}, 'linewidth', 2);
            end
            
            % CSF fitting with smoothing spline (Method 2).
            if (FitSmoothSpline)
                if (OneFigurePerSub)
                    plot(smoothPlotSFVals,smoothPlotPreds,colorOptionsCSF2{ff},'LineWidth',4);
                else
                    plot(smoothPlotSFVals,smoothPlotPreds,'c-','LineWidth',4);
                end
            end
            
            % Add details per each plot of the subject.
            if (~OneFigurePerSub)
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
            if (~OneFigurePerSub)
                if (WaitForKeyToPlot)
                    disp('Press a key to draw next plot!');
                    pause;
                    close all;
                end
            end
        end
        
        %% Add details per each plot of the subject.
        if (OneFigurePerSub)
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
            if (WaitForKeyToPlot)
                disp('Press a key to draw next plot!');
                pause;
                close all;
            end
        end
    end
end
