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
%    02/13/23   smo    - Now we fit and plot the data with all methods at
%                        the same time.
%    02/15/23   smo    - Added an option to save the CSF plot.
%    02/21/23   smo    - Added an option to calculate and bootstrap AUC.
%                        Also, we can choose which domain to fit CSF either
%                        linear or log.

%% Initialize.
clear; close all;

%% Decide the options to fit CSF and how to plot.
%
% Plotting options.
OneFigurePerSub = false;
WaitForKeyToPlot = true;
SaveCSFPlot = false;
PlotAUC = true;

% Fitting options.
FitAsymmetricParabolic = false;
FitSmoothSpline = true;
CalAUC = true;
BootstrapAUC = false;
CSFFittingDomain = 'log';

% When fitting Smooth spline, You can choose option among {'crossVal',
% 'crossValBootWithin', 'crossValBootAcross', 'type'}.
if (FitSmoothSpline)
    optionSearchSmoothParamSet = {'crossValBootAcross'};
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
            figurePosition = [200 300 figureSize figureSize];
            set(gcf,'position',figurePosition);
            
            % Cross-validation figure info.
            if (FitSmoothSpline)
                crossFig = figure; hold on;
                figurePositionCross = [200+figureSize 300 figureSize figureSize];
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
                figurePosition = [200 300 figureSize figureSize];
                set(gcf,'position',figurePosition);
                
                % Cross-validation figure info.
                if (FitSmoothSpline)
                    crossFig = figure; hold on;
                    figurePositionCross = [200+figureSize 300 figureSize figureSize];
                    set(gcf,'position',figurePositionCross);
                end
            end
            
            % Read out the variables per each filter. These values are
            % linear units, which we will conver them on log space.
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
            
            %% Calculate log sensitivity.
            %
            % Raw data.
            sensitivity = log10(1./thresholds);
            sensitivityMedianBoot = log10(1./medianThresholdsBoot);
            
            % Bootstrapped values.
            %
            % For calculation of confindence
            % interval from bootstrap, (low) threshold becomes (high)
            % sensitivity, and vice versa.
            sensitivityBoot = log10(1./squeeze(thresholdsBoot));
            sensitivityBootHigh = log10(1./lowThresholdBoot);
            sensitivityBootLow = log10(1./highThresholdBoot);
            
            % Additional bootstrapped values for cross-validation.
            sensitivityBootCross1 = log10(1./squeeze(thresholdsBootCross1));
            sensitivityBootCross2 = log10(1./squeeze(thresholdsBootCross2));
            
            % Calculate spatial frequency in log space.
            sineFreqCyclesPerDegLog = log10(sineFreqCyclesPerDegNum);
            
            %% Sort each array in an ascending order of spatial frequency.
            [sineFreqCyclesPerDegLogSorted I] = sort(sineFreqCyclesPerDegLog,'ascend');
            
            % Sorted according to the order of spatial frequency.
            sensitivitySorted = sensitivity(I);
            sensitivityMedianBootSorted = sensitivityMedianBoot(I);
            sensitivityBootHighSorted = sensitivityBootHigh(I);
            sensitivityBootLowSorted = sensitivityBootLow(I);
            sensitivityBootSorted = sensitivityBoot(I,:);
            sensitivityBootCross1Sorted = sensitivityBootCross1(I,:);
            sensitivityBootCross2Sorted = sensitivityBootCross2(I,:);
            
            %% Set variables to fit CSF.
            switch CSFFittingDomain
                case 'log'
                    mySFVals = sineFreqCyclesPerDegLogSorted;
                    myCSVals = sensitivitySorted;
                    myWs = 1./(sensitivityBootHighSorted-sensitivityBootLowSorted);
                case 'linear'
                    mySFVals = 10.^sineFreqCyclesPerDegLogSorted;
                    myCSVals = 10.^sensitivitySorted;
                    myWs = 1./10.^(sensitivityBootHighSorted-sensitivityBootLowSorted);
            end
            
            %% Fitting method 1) Asymmetric parabolic function.
            if (FitAsymmetricParabolic)
                % Set parameters for optimazation of the parameter p.
                p0 = [log10(10) log10(0.5) 0.1 0.1];
                p_lowerBound = [log10(10) log10(0.5) 0.1 0.1];
                p_higherBound = [log10(400) log10(18) 10 10];
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
                myCSValsBoot = sensitivityBootSorted';
                myCSValsCross1 = sensitivityBootCross1Sorted';
                myCSValsCross2 = sensitivityBootCross2Sorted';
                nBootPoints = size(myCSValsBoot,1);
                
                %% Search smoothing parameter here.
                %
                % Set the number of points for plotting the results.
                nSmoothPoints = 100;
                
                % Set the smoothing param searching options.
                nSmoothingParams = 30;
                minSmoothingParam = 0;
                maxSmoothingParam = 1;
                crossSmoothingParams = linspace(minSmoothingParam,maxSmoothingParam,nSmoothingParams);
                
                % We search smoothing parameter based on the option selected
                % from the abovOptionSearchSmoothParamSete.
                nOptionsSearchSmoothParamSet = length(optionSearchSmoothParamSet);
                
                % Make a loop here to test all the methods selected from
                % the above.
                for oo = 1:nOptionsSearchSmoothParamSet
                    OptionSearchSmoothParam = optionSearchSmoothParamSet{oo};
                    
                    % Switch the method how to search smoothing parameter.
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
                                if (sss == round(length(crossSmoothingParams)*0.25))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam, '25%');
                                elseif (sss == round(length(crossSmoothingParams)*0.50))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '50%');
                                elseif (sss == round(length(crossSmoothingParams)*0.75))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam, '75%');
                                elseif (sss == length(crossSmoothingParams))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam, '100%');
                                end
                            end
                            
                        case 'crossValBootAcross'
                            % Make a loop for bootstrapping AUC if you
                            % want. This should start before deciding the
                            % searching method for smoothing param, but for
                            % now we start from here as we are leaning to
                            % use this method.
                            if (BootstrapAUC)
                                nBootstrapAUC = 10;
                            else
                                nBootstrapAUC = 1;
                            end
                            
                            for aaa = 1:nBootstrapAUC
                                
                                % Make a loop for testing smoothing paramemters.
                                for sss = 1:length(crossSmoothingParams)
                                    smoothCrossError(sss,aaa) = 0;
                                    
                                    % Bootstrap for cross-validation happens here.
                                    nCrossValBootAcross = 20;
                                    for cc = 1:nCrossValBootAcross
                                        
                                        % Draw new fit/cross dataset (N=20)
                                        % out of bootstrapped values
                                        % (N=100). Once we drew the set, we
                                        % will use the same set for all
                                        % smoothing params.
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
                                        smoothCrossError(sss,aaa) = smoothCrossError(sss,aaa) + sum((bootCSFDataCross{cc}' - smoothDataPredsCross).^2);
                                    end
                                    
                                    % Print out the progress.
                                    if (sss == round(length(crossSmoothingParams)*0.25))
                                        fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam, '25%');
                                    elseif (sss == round(length(crossSmoothingParams)*0.50))
                                        fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '50%');
                                    elseif (sss == round(length(crossSmoothingParams)*0.75))
                                        fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '75%');
                                    elseif (sss == length(crossSmoothingParams))
                                        fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '100%');
                                    end
                                end
                                % Print out the progress of bootstrapping
                                % AUC. It will take a while.
                                fprintf('\t Bootstrapping AUC progress - (%d/%d) \n', aaa, nBootstrapAUC);
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
                                
                                % Print out the progress.
                                if (sss == round(length(crossSmoothingParams)*0.25))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam, '25%');
                                elseif (sss == round(length(crossSmoothingParams)*0.50))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '50%');
                                elseif (sss == round(length(crossSmoothingParams)*0.75))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '75%');
                                elseif (sss == length(crossSmoothingParams))
                                    fprintf('Method = (%s) / Smoothing param testing progress - (%s) \n', OptionSearchSmoothParam,  '100%');
                                end
                            end
                            
                        case 'type'
                            % Type a number manually.
                            smoothingParam = 0.1;
                    end
                    
                    %% Plot the bootstrapped results if we did.
                    %
                    % Make a loop to repeat this part to bootstrap AUC.
                    for aaa = 1:nBootstrapAUC
                        
                        if any(strcmp(OptionSearchSmoothParam,{'crossValBootWithin','crossValBootAcross','crossVal'}))
                            % Set the smoothing params that has the smallest error.
                            [~,index] = min(smoothCrossError(:,aaa));
                            smoothingParam = crossSmoothingParams(index);
                            
                            % Plot cross validation error.
                            %
                            % If we test multiple options at the same time, we
                            % only make new figure at the beginning.
                            if (nOptionsSearchSmoothParamSet > 1)
                                if (oo == 1)
                                    figure(crossFig); hold on;
                                end
                            else
                                figure(crossFig); hold on;
                            end
                            
                            % If we make a loop for bootstrapping AUC, we
                            % plot this part only for the very first set.
                            % This is the same for plotting the CSF.
                            if (aaa == 1)
                                markerColorOptionsSmoothSpline = {'r','g','b','c'};
                                plot(crossSmoothingParams,smoothCrossError(:,aaa),'ko','MarkerSize',6);
                                plot(smoothingParam,smoothCrossError(index,aaa),'co','MarkerSize',8,'Markerfacecolor',markerColorOptionsSmoothSpline{oo},'Markeredgecolor','k');
                                xlabel('Smoothing parameter','fontsize',15);
                                ylabel('Cross-validation errors','fontsize',15);
                                title('Cross-validation error accoring to smoothing parameter','fontsize',15);
                                xlim([minSmoothingParam maxSmoothingParam]);
                            end
                        end
                        
                        % Get the values for plotting smooth spline fitting curve.
                        smoothFit = fit(mySFVals',myCSVals','smoothingspline','SmoothingParam',smoothingParam);
                        smoothPlotSFVals{oo,aaa} = log10(logspace(min(mySFVals),max(mySFVals),nSmoothPoints))';
                        smoothPlotPreds{oo,aaa} = feval(smoothFit,smoothPlotSFVals{oo,aaa});
                        
                        %% Get the area under the CSF curve (AUC).
                        if (CalAUC)
                            nPointsCalAUC = 1000;
                            calAUCSFVals{oo} = log10(logspace(min(mySFVals),max(mySFVals),nPointsCalAUC))';
                            calAUCPreds{oo} = feval(smoothFit,calAUCSFVals{oo});
                            AUC(aaa) = 0;
                            for aa = 1:nPointsCalAUC-1
                                AUCTemp = (calAUCSFVals{oo}(aa+1)-calAUCSFVals{oo}(aa)) * calAUCPreds{oo}(aa);
                                AUC(aaa) = AUC(aaa) + AUCTemp;
                            end
                            fprintf('Calculated AUC (%d/%d) is (%.5f) \n',...
                                aaa, nBootstrapAUC, AUC(aaa));
                            meanAUC = mean(AUC);
                            stdAUC = std(AUC);
                        end
                    end
                end
                
                % Add legend to cross figure.
                f_cross = flip(get(gca, 'Children'));
                idxLegendCross = [2:2:nOptionsSearchSmoothParamSet*2];
                legend(f_cross(idxLegendCross),optionSearchSmoothParamSet,'fontsize',12);
            end
            
            %% Plot data figure here.
            figure(dataFig);
            
            % Set marker/line options for the plot.
            colorOptionsRaw = {'k.','r.','g.','b.','c.'};
            colorOptionsBoot = {'k+','r+','g+','b+','c+'};
            colorOptionsCSF = {'k-','r-','g-','b-','c-'};
            colorOptionsCSF2 = {'k--','r--','g--','b--','c--'};
            colorOptionsCI = {'k','r','g','b','c'};
            colorOptionsSmoothSpline = {'r-','g--','b:'};
            
            % Raw data.
            switch CSFFittingDomain
                case 'log'
                    plot(sineFreqCyclesPerDegLogSorted, sensitivitySorted, colorOptionsRaw{ff},'markersize',20);
                    plot(sineFreqCyclesPerDegLogSorted, sensitivityMedianBootSorted, colorOptionsBoot{ff},'markersize',20);
                case 'linear'
                    plot(10.^sineFreqCyclesPerDegLogSorted, 10.^sensitivitySorted, colorOptionsRaw{ff},'markersize',20);
                    plot(10.^sineFreqCyclesPerDegLogSorted, 10.^sensitivityMedianBootSorted, colorOptionsBoot{ff},'markersize',20);
            end
            
            % Confidence Interval.
            switch CSFFittingDomain
                case 'log'
                    errorNeg = abs(sensitivityMedianBootSorted - sensitivityBootLowSorted);
                    errorPos = abs(sensitivityBootHighSorted - sensitivityMedianBootSorted);
                    e = errorbar(sineFreqCyclesPerDegLogSorted, sensitivityMedianBootSorted, ...
                        errorNeg, errorPos, colorOptionsCI{ff});
                case 'linear'
                    errorNeg = abs(10.^sensitivityMedianBootSorted - 10.^sensitivityBootLowSorted);
                    errorPos = abs(10.^sensitivityBootHighSorted - 10.^sensitivityMedianBootSorted);
                    e = errorbar(10.^sineFreqCyclesPerDegLogSorted, 10.^sensitivityMedianBootSorted, ...
                        errorNeg, errorPos, colorOptionsCI{ff});
            end
            e.LineStyle = 'none';
            
            % CSF fitting results with asymmetric parabolic (Method 1).
            if (FitAsymmetricParabolic)
                SF_CSF_start = log10(3);
                SF_CSF_end = log10(18);
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
                    for oo = 1:nOptionsSearchSmoothParamSet
                        % We will plot the CSF derived from the very first
                        % set if we do bootstrapping for AUC.
                        numBootstrapAUCToPlotCSF = 1;
                        plot(smoothPlotSFVals{oo,numBootstrapAUCToPlotCSF},smoothPlotPreds{oo,numBootstrapAUCToPlotCSF},...
                            colorOptionsSmoothSpline{oo},'LineWidth',4);
                    end
                end
            end
            
            % Plot AUC results if you want.
            if (PlotAUC)
                for aa = 1:nPointsCalAUC
                    plot(ones(1,2)*calAUCSFVals{1}(aa), [0 calAUCPreds{1}(aa)],'color',[1 0 0 0.1]);
                end
            end
            
            % Add details per each plot of the subject.
            if (~OneFigurePerSub)
                xlabel('Spatial Frequency (cpd)','fontsize',15);
                ylabel('Contrast Sensitivity','fontsize',15);
                
                switch CSFFittingDomain
                    case 'log'
                        xticks(sineFreqCyclesPerDegLogSorted);
                        xticklabels(10.^sineFreqCyclesPerDegLogSorted);
                        yaxisRange = log10([0:50:300]);
                        ylim(log10([1 300]));
                        yticks(yaxisRange);
                        yticklabels(10.^yaxisRange);
                    case 'linear'
                        xticks(10.^sineFreqCyclesPerDegLogSorted);
                        xticklabels(10.^sineFreqCyclesPerDegLogSorted);
                        yaxisRange = ([0:50:300]);
                        ylim([0 300]);
                        yticks(yaxisRange);
                        yticklabels(yaxisRange);
                end
                
                title(sprintf('CSF curve - Sub %s',subjectName),'fontsize',15);
                subtitle('Fitting was done on log-log space');
                
                % Add legend.
                f_data = flip(get(gca, 'Children'));
                
                numSpaceLegend = 4;
                idxLegendRaw = linspace(1, 1+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
                idxLegendBoot = linspace(2, 2+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
                
                filterWls = {'392 nm' '417 nm' '437 nm' '456 nm' '476 nm'};
                for ll = 1:nSineFreqCyclesPerDeg
                    contentLegendRaw{ll} = sprintf('%s (%s) -PF',theData.filterOptions{ll}, filterWls{ll});
                    contentLegendBoot{ll} = sprintf('%s (%s) -Boot',theData.filterOptions{ll}, filterWls{ll});
                end
                
                % Add legend when drawing one figure per each filter.
                legend(f_data([1,2,4:end]),[append(filterWls{ff},'-PF'), append(filterWls{ff},'-Boot'), ...
                    optionSearchSmoothParamSet],'fontsize',13,'location', 'northeast');
                
                % Add AUC to the plot.
                textAUC = sprintf('AUC = %.4f (%.4f)', meanAUC, stdAUC);
                switch CSFFittingDomain
                    case 'log'
                        text(log10(3),log10(1.5),textAUC,'color','r','fontsize',15);
                    case 'linear'
                        text(3,5,textAUC,'color','r','fontsize',15);
                end
            end
            
            % Save the CSF plot if you want.
            if (SaveCSFPlot)
                if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),...
                        subjectName,'CSF');
                    filterNames = {'A', 'B', 'C', 'D', 'E'};
                    testFilename = fullfile(testFiledir, sprintf('%s_%s_%s','CSF', subjectName, filterNames{ff}));
                    testFileFormat = '.tiff';
                    saveas(dataFig, append(testFilename,testFileFormat));
                    disp('CSF plot has been saved successfully!');
                end
            end
            
            % Key stroke to start measurement.
            if (~OneFigurePerSub)
                if (WaitForKeyToPlot)
                    fprintf('\t Press a key to draw next plot! \n');
                    pause;
                    close all;
                end
            end
        end
        
        %% Add details per each plot of the subject.
        if (OneFigurePerSub)
            xlabel('Spatial Frequency (cpd)','fontsize',15);
            ylabel('Contrast Sensitivity','fontsize',15);
            
            xticks(sineFreqCyclesPerDegLogSorted);
            xticklabels(sineFreqCyclesPerDegLogSorted);
            
            yaxisRange = [0:50:300];
            ylim([0 300]);
            yticks(yaxisRange);
            yticklabels(yaxisRange);
            
            title(sprintf('CSF curve - Sub %s',subjectName),'fontsize',15);
            
            % Add legend.
            f_data = flip(get(gca, 'Children'));
            
            numSpaceLegend = 4;
            idxLegendRaw = linspace(1, 1+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            idxLegendBoot = linspace(2, 2+numSpaceLegend*(nSineFreqCyclesPerDeg-1), nSineFreqCyclesPerDeg);
            
            filterWls = {'392 nm' '417 nm' '437 nm' '456 nm' '476 nm'};
            for ll = 1:nSineFreqCyclesPerDeg
                contentLegendRaw{ll} = sprintf('%s (%s) -PF',theData.filterOptions{ll}, filterWls{ll});
                contentLegendBoot{ll} = sprintf('%s (%s) -Boot',theData.filterOptions{ll}, filterWls{ll});
            end
            
            % Add legend when drawing one figure per each subject.
            legend(f_data([idxLegendRaw idxLegendBoot]), [contentLegendRaw contentLegendBoot], ...
                'fontsize', 13, 'location', 'northeast');
            
            % Save the CSF plot if you want.
            if (SaveCSFPlot)
                if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),...
                        subjectName,'CSF');
                    filterNames = {'A', 'B', 'C', 'D', 'E'};
                    testFilename = fullfile(testFiledir, sprintf('%s_%s_%s','CSF', subjectName, filterNames{ff}));
                    testFileFormat = '.tiff';
                    saveas(dataFig, append(testFilename,testFileFormat));
                    disp('CSF plot has been saved successfully!');
                end
            end
            
            % Key stroke to start measurement.
            if (WaitForKeyToPlot)
                disp('Press a key to draw next plot!');
                pause;
                close all;
            end
        end
    end
end
