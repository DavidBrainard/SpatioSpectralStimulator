function [paramsFitted, ...
    thresholdFitted, thresholdFittedBoot, medianThresholdBoot, lowThresholdBoot, highThresholdBoot,...
    slopeFitted,medianSlopeBoot,lowSlopeBoot,highSlopeBoot,thresholdFittedBootCross1,thresholdFittedBootCross2, ...
    legendHandles] = FitPFToData(stimLevels,pCorrect,options)
% Fit Psychometric function to the given data.
%
% Syntax:
%    [paramsFitted,...
%      thresholdFitted, medianThresholdBoot,lowThresholdBoot,highSlopeBoot ...
%      slopeFitted,medianSlopeBoot,lowSlopeBoot,highSlopeBoot] = FitPFToData(stimLevels,pCorrect)
%
% Description:
%    This fits Psychometric function to the data given as an input. You can
%    choose PF either Weibull or Logistic.
%
% Inputs:
%    stimLevels -                 Array of the stimulus levels.
%    pCorrect -                   Array of percentage correct per each
%                                 stimulus level. This should be the same
%                                 size of the stimLevels.
%
% Outputs:
%    paramsFitted -               Parameters found from PF fitting. It's in
%                                 the format of the array [threshold slope
%                                 guess lapse]. You can choose which one to
%                                 be free/not free parameters.
%    thresholdFitted -            Threshold at criterion
%    thresholdFittedBoot -        All bootstrapped values.
%    medianThresholdBoot -        Median of boostrapped threshold
%    lowThresholdBoot -           Low end of bootstrapped CI
%    highThresholdBoot -          High end of bootstrapped CI
%    slopeFitted -                Fit slope
%    medianSlopeBoot -            Median of bootstrapped slopes
%    lowSlopeBoot -               Low end of bootstrapped slope CI
%    highSlopeBoot -              High end of bootstrapped slope CI
%    thresholdFittedBootCross1 -  Additional bootstrapped values for
%                                 cross-validation - training.
%    thresholdFittedBootCross2 -  Additional bootstrapped values for
%                                 cross-validation - validating.
%    legendHandles -              Vector of plot handles so calling routine
%                                 can make a sensible legend.
%
% Optional key/value pairs:
%    PF -                         Default to 'weibull'. Choose the function
%                                 to fit the data either 'weibull' or
%                                 'logistic'.
%    paramsFree -                 Default to [1 1 0 1]. This decides which
%                                 parameters to be free. Array represents
%                                 [threshold slope guess lapse]. Each can
%                                 be set either 0 or 1, where 0 = fixed
%                                 1 = free.
%    beta -                       Set of slope parameters to fit for.  Each
%                                 of these is fixed and the best fit over
%                                 the set is used.  This allows us to
%                                 control the slope range. When this is not
%                                 empty, the second entry of params free
%                                 above is ignored, and the fit is done for
%                                 the passed set with best returned.
%                                 Default is empty.
%    nTrials -                    Default to 20. The number of trials
%                                 conducted per each stimulus level.
%    thresholdCriterion -         Default to 0.81606. This is the value of
%                                 pCorrect as a criteria to find threshold.
%    newFigureWindow -            Default to true. Make a new figure window
%                                 to plot the results. If you want to plot
%                                 multiple threshold results in a subplot,
%                                 set this to false and do so.
%    pointSize                  - Default to 100 each. Set the size of data
%                                 point on the scatter plot.
%    axisLog                    - Default to true. If it sets to true, plot
%                                 the graph with x-axis on log space.
%    questPara                  - Default to blank. Add Quest fit if this
%                                 is not empty.Row vector or matrix of
%                                 parameters.
%                                 threshold  Threshold in log unit
%                                 slope      Slope
%                                 guess      Guess rate
%                                 lapse      Lapse rate
%                                 Parameterization matches the Mathematica
%                                 code from the Watson QUEST+ paper.
%    addLegend                  - Default to true. Add legend when it sets
%                                 to true.
%    nBootstraps                - Number of bootstraps to run.  0 means no
%                                 bootstrapping. Default 0.
%    bootConfInterval           - Size of bootstrapped confidence interval.
%                                 Default 0.8.
%    verbose -                  - Default to true. Boolean. Controls
%                                 plotting and printout.
%
% See also:
%    N/A

% History:
%   02/25/22  dhb, smo         - Started on it.
%   03/14/22  smo              - Added a plotting option to make different
%                                marker size over different number of
%                                trials.
%   11/14/22  dhb              - Bootstrapping
%   02/06/23  smo              - Now print out all bootstrap values.
%   04/03/23  smo              - Modified to prevent generating negative
%                                threshold values on linear space when
%                                bootstrapping.
%   04/10/23  smo              - Deleted the while loop and set the lowest
%                                possible threshold to 0.00001.

%% Set parameters.
arguments
    stimLevels
    pCorrect
    options.PF = 'weibull'
    options.paramsFree (1,4) = [1 1 0 1]
    options.beta = []
    options.nTrials (1,1) = 20
    options.thresholdCriterion (1,1) = 0.81606
    options.newFigureWindow (1,1) = true
    options.pointSize = ones(1,length(stimLevels))*100
    options.axisLog (1,1) = true
    options.questPara = []
    options.addLegend (1,1) = true
    options.nBootstraps = 0
    options.bootConfInterval = 0.8
    options.verbose (1,1) = true
end

%% Check the size of the input parameters.
if (~any(size(stimLevels) == size(pCorrect)))
    error('Stimulus level and pCorrect array size does not match!');
end

%% Set PF fitting type here.
switch options.PF
    case 'weibull'
        PF = @PAL_Weibull;
    case 'logistic'
        PF = @PAL_Logistic;
    otherwise
end

%% Fitting PF.
%
% Set up the PF fitting (requires Palamedes toolbox).  Note that the
% catch trials are added in here in the call to the fit.
nTrialsPerContrast = options.nTrials * ones(size(stimLevels));
nCorrect = round(pCorrect .* nTrialsPerContrast);

% Set an initial search parameters, with
% gridded slope (aka beta).  The grid search
% is only done if options.beta is empty.
searchGrid.alpha = mean(stimLevels);
searchGrid.beta = 10.^(-2:0.01:2);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.01;
lapseLimits = [0 0.05];

% PF fitting happens here. Search over passed
% list of possible slopes
if (~isempty(options.beta))
    paramsFreeUse = options.paramsFree;
    paramsFreeUse(2) = 0;
    for ss = 1:length(options.beta)
        searchGridUse = searchGrid;
        searchGridUse.beta = options.beta(ss);
        [paramsFittedList(ss,:) LL(ss)] = PAL_PFML_Fit(stimLevels, nCorrect, ...
            nTrialsPerContrast, searchGridUse, paramsFreeUse, PF, 'lapseLimits', lapseLimits);
    end
    [~,index] = max(LL);
    paramsFitted = paramsFittedList(index,:);
    PLOT_SLOPELL = false;
    if (PLOT_SLOPELL)
        theFig = gcf;
        llfig = figure; clf; hold on
        plot(options.beta,LL,'ro');
        figure(theFig);
        pause;
        close(llfig);
    end
    
    % Otherwise use Palamedes grid search
else
    paramsFitted = PAL_PFML_Fit(stimLevels, nCorrect, ...
        nTrialsPerContrast, searchGrid, options.paramsFree, PF, 'lapseLimits', lapseLimits);
end
thresholdFitted = PF(paramsFitted, options.thresholdCriterion, 'inv');
slopeFitted = paramsFitted(2);

%% Bootstrap fits.
if (options.nBootstraps > 0)
    paramsFittedBoot = zeros(options.nBootstraps,4);
    for bb = 1:options.nBootstraps
        % Set the array size for bootstrapped data.
        zeroBaseArrayBoot = zeros(size(nCorrect));
        
        % Set an array for saving the bootstrapped data.
        nCorrectBoot = zeroBaseArrayBoot;
        
        % Additional bootstrap data arrays for cross-validation for fitting
        % CSF. We make two more.
        nCorrectCross1 = zeroBaseArrayBoot;
        nCorrectCross2 = zeroBaseArrayBoot;
        nTrialsPerContrastCross1 = zeroBaseArrayBoot;
        nTrialsPerContrastCross2 = zeroBaseArrayBoot;
        
        for cc = 1:length(nTrialsPerContrast)
            % Generate the bootstrapped data.
            trialsBoot = zeros(1,nTrialsPerContrast(cc));
            trialsBoot(1:nCorrect(cc)) = 1;
            index = randi(nTrialsPerContrast(cc),1,nTrialsPerContrast(cc));
            nCorrectBoot(cc) = sum(trialsBoot(index));
            
            % Here we shuffle the bootstrapped data to generate the
            % additional data. Here we make two more separate dataset, each
            % uses the half size of trials per each contrast.
            %
            % Therefore, if the original dataset uses a total of 20 trials
            % per each contrast, these additional data will use 10 trials.
            trialsShuffle = Shuffle(trialsBoot);
            splitN = round(length(trialsShuffle)/2);
            
            nCorrectCross1(cc) = sum(trialsShuffle(1:splitN));
            nCorrectCross2(cc) = sum(trialsShuffle((splitN+1):end));
            nTrialsPerContrastCross1(cc) = length(trialsShuffle(1:splitN));
            nTrialsPerContrastCross2(cc) = length(trialsShuffle((splitN+1):end));
        end
        
        %% Fit the bootstrap data in the same way we fit actual data.
        %
        % 1) Main bootstrapped data.
        if (~isempty(options.beta))
            paramsFreeUse = options.paramsFree;
            paramsFreeUse(2) = 0;
            for ss = 1:length(options.beta)
                searchGridUse = searchGrid;
                searchGridUse.beta = options.beta(ss);
                [paramsFittedList(ss,:) LL(ss)] = PAL_PFML_Fit(stimLevels, nCorrectBoot, ...
                    nTrialsPerContrast, searchGridUse, paramsFreeUse, PF, 'lapseLimits', lapseLimits);
            end
            [~,index] = max(LL);
            paramsFittedBoot(bb,:) = paramsFittedList(index,:);
        else
            paramsFittedBoot(bb,:) = PAL_PFML_Fit(stimLevels, nCorrectBoot, ...
                nTrialsPerContrast, searchGrid, options.paramsFree, PF, 'lapseLimits', lapseLimits);
        end
        
        % Grab bootstrapped threshold.
        thresholdFittedBoot(bb) = PF(paramsFittedBoot(bb,:), options.thresholdCriterion, 'inv');
        
        % Limit the lowest contrast threshold value here. For some cases,
        % it gives the negative threshold results and we will convert them
        lowLimitContrastThresholdLinear = 0.00001;
        thresholdFittedBoot(thresholdFittedBoot<lowLimitContrastThresholdLinear) = lowLimitContrastThresholdLinear;
        
        % 2) Cross data 1.
        if (~isempty(options.beta))
            paramsFreeUse = options.paramsFree;
            paramsFreeUse(2) = 0;
            for ss = 1:length(options.beta)
                searchGridUse = searchGrid;
                searchGridUse.beta = options.beta(ss);
                [paramsFittedList(ss,:) LL(ss)] = PAL_PFML_Fit(stimLevels, nCorrectCross1, ...
                    nTrialsPerContrastCross1, searchGridUse, paramsFreeUse, PF, 'lapseLimits', lapseLimits);
            end
            [~,index] = max(LL);
            paramsFittedBootCross1(bb,:) = paramsFittedList(index,:);
        else
            paramsFittedBootCross1(bb,:) = PAL_PFML_Fit(stimLevels, nCorrectCross1, ...
                nTrialsPerContrastCross1, searchGrid, options.paramsFree, PF, 'lapseLimits', lapseLimits);
        end
        
        % Grab bootstrapped threshold.
        thresholdFittedBootCross1(bb) = PF(paramsFittedBootCross1(bb,:), options.thresholdCriterion, 'inv');
        
        % Limit the lowest contrast threshold value.
        thresholdFittedBootCross1(thresholdFittedBootCross1<lowLimitContrastThresholdLinear) = lowLimitContrastThresholdLinear;
        
        % 3) Cross data 2.
        if (~isempty(options.beta))
            paramsFreeUse = options.paramsFree;
            paramsFreeUse(2) = 0;
            for ss = 1:length(options.beta)
                searchGridUse = searchGrid;
                searchGridUse.beta = options.beta(ss);
                [paramsFittedList(ss,:) LL(ss)] = PAL_PFML_Fit(stimLevels, nCorrectCross2, ...
                    nTrialsPerContrastCross2, searchGridUse, paramsFreeUse, PF, 'lapseLimits', lapseLimits);
            end
            [~,index] = max(LL);
            paramsFittedBootCross2(bb,:) = paramsFittedList(index,:);
        else
            paramsFittedBootCross2(bb,:) = PAL_PFML_Fit(stimLevels, nCorrectCross2, ...
                nTrialsPerContrastCross2, searchGrid, options.paramsFree, PF, 'lapseLimits', lapseLimits);
        end
        
        % Grab bootstrapped threshold.
        thresholdFittedBootCross2(bb) = PF(paramsFittedBootCross2(bb,:), options.thresholdCriterion, 'inv');
        
        % Limit the lowest contrast threshold value.
        thresholdFittedBootCross2(thresholdFittedBootCross2<lowLimitContrastThresholdLinear) = lowLimitContrastThresholdLinear;
    end
    
    % Extract some useful data from main bootstrapped data.
    medianThresholdBoot = median(thresholdFittedBoot);
    lowThresholdBoot = prctile(thresholdFittedBoot,100*(1-options.bootConfInterval)/2);
    highThresholdBoot = prctile(thresholdFittedBoot,100-100*(1-options.bootConfInterval)/2);
    medianSlopeBoot = median(paramsFittedBoot(:,2));
    lowSlopeBoot = prctile(paramsFittedBoot(:,2),100*(1-options.bootConfInterval)/2);
    highSlopeBoot = prctile(paramsFittedBoot(:,2),100-100*(1-options.bootConfInterval)/2);
    
else
    paramsFittedBoot = [];
    thresholdFittedBoot = [];
    medianThresholdBoot = [];
    lowThresholdBoot = [];
    highThresholdBoot = [];
    medianSlopeBoot = [];
    lowSlopeBoot = [];
    highSlopeBoot = [];
    thresholdFittedBootCross1 = [];
    thresholdFittedBootCross2 = [];
end

% Make a smooth curves with finer stimulus levels.
nFineStimLevels = 1000;
fineStimLevels = linspace(0, max(stimLevels), nFineStimLevels);
smoothPsychometric = PF(paramsFitted, fineStimLevels);

if (options.verbose)
    fprintf('Threshold was found at %.4f (linear unit) \n', thresholdFitted);
end

%% Plot the results if you want.
if (options.verbose)
    if (options.newFigureWindow)
        figure; clf; hold on;
    end
    
    % Plot all experimental data (gray points).
    %
    % Marker size will be different over the number of the trials per each
    % test point.
    %
    % Plot it on log space if you want.
    if (options.axisLog)
        stimLevelsPlot = log10(stimLevels);
        thresholdFittedLog = log10(thresholdFitted);
        fineStimLevelsPlot = log10(fineStimLevels);
    else
        stimLevelsPlot = stimLevels;
        fineStimLevelsPlot = fineStimLevels;
    end
    
    % Plot bootstraps
    for bb = 1:options.nBootstraps
        % We will skip plotting if threshold was found at negative values.
        if (paramsFittedBoot(bb,1)<0)
            continue;
        end
        smoothPsychometricBoot = PF(paramsFittedBoot(bb,:), fineStimLevels);
        h_bsfit = plot(fineStimLevelsPlot,smoothPsychometricBoot,'Color',[0.9 0.8 0.8],'LineWidth',0.5);
    end
    
    % Plot best fit here.
    %
    % Raw data.
    h_data = scatter(stimLevelsPlot, pCorrect, options.pointSize,...
        'MarkerEdgeColor', zeros(1,3), 'MarkerFaceColor', ones(1,3) * 0.5, 'MarkerFaceAlpha', 0.5);
    
    % Plot the PF fitting results.
    %
    % Set the line color differently if it passes lower than expecteced
    % number of fitting points. We expect to have the number of 8 points to
    % fit for SACC project, but as we re-fit the PF without bad contrast
    % points, so we want to plot it in different color. We may delete this
    % part later on.
    %
    % Now it always draw in red color (as of 10/12/23).
    PFLineColor = 'r';
    h_pffit = plot(fineStimLevelsPlot,smoothPsychometric,PFLineColor,'LineWidth',3);
    
    % Mark the threshold point (red point or yellow).
    if(options.axisLog)
        threshMarkerFaceColor = 'r';
        h_thresh = plot(thresholdFittedLog,options.thresholdCriterion,'ko','MarkerFaceColor',threshMarkerFaceColor,'MarkerSize',12);
        
        % Mark the median bootstrapped threshold point (green point) and
        % its confidence interval (green line).
        if (options.nBootstraps > 0)
            h_bsthresh = errorbarX(log10(medianThresholdBoot),options.thresholdCriterion,log10(medianThresholdBoot)-log10(lowThresholdBoot),log10(highThresholdBoot)-log10(medianThresholdBoot),'go');
            set(h_bsthresh,'MarkerSize',9); set(h_bsthresh,'MarkerFaceColor',[0.1 0.7 0.1]); set(h_bsthresh,'MarkerEdgeColor','k'); set(h_bsthresh(2),'LineWidth',0.5); set(h_bsthresh(1),'LineWidth',3); set(h_bsthresh(1),'color',[0.1 0.7 0.1]);
        end
        xlabel('Log Contrast', 'FontSize', 15);
    else
        h_thresh = plot(thresholdFitted,options.thresholdCriterion,'ko','MarkerFaceColor','r','MarkerSize',12);
        if (options.nBootstraps > 0)
            h_bsthresh = errorbarX(medianThresholdBoot,options.thresholdCriterion,medianThresholdBoot-lowThresholdBoot,highThresholdBoot-medianThresholdBoot,'go');
            set(h_bsthresh,'MarkerSize',9); set(h_bsthresh,'MarkerFaceColor',[0.1 0.7 0.1]); set(h_bsthresh,'MarkerEdgeColor','k'); set(h_bsthresh(2),'LineWidth',0.5); set(h_bsthresh(1),'LineWidth',3); set(h_bsthresh(1),'color',[0.1 0.7 0.1]);
        end
        xlabel('Contrast', 'FontSize', 15);
    end
    ylabel('pCorrect', 'FontSize', 15);
    ylim([0 1]);
    if (~isempty(options.beta))
        if (options.nBootstraps > 0)
            title({sprintf('Threshold: %0.4f, [%0.4f (%0.4f, %0.4f, %d%%]',...
                thresholdFitted,medianThresholdBoot,lowThresholdBoot,highThresholdBoot,round(100*options.bootConfInterval)) ; ...
                sprintf('Slope range: %0.4f to %0.4f',min(options.beta),max(options.beta)); ...
                sprintf('Slope: %0.4f, [%0.4f (%0.4f, %0.4f, %d%%]',...
                slopeFitted,medianSlopeBoot,lowSlopeBoot,highSlopeBoot,round(100*options.bootConfInterval))});
        else
            title({sprintf('Threshold: %0.4f',thresholdFitted); sprintf('Slope: %0.4f',slopeFitted); ...
                sprintf('Slope range: %0.4f to %0.4f',min(options.beta),max(options.beta))})
        end
    else
        if (options.nBootstraps > 0)
            title({sprintf('Threshold: %0.4f, [%0.4f (%0.4f, %0.4f, %d%%]',...
                thresholdFitted,medianThresholdBoot,lowThresholdBoot,highThresholdBoot,round(100*options.bootConfInterval)) ; ...
                sprintf('Slope: %0.4f, [%0.4f (%0.4f, %0.4f, %d%%]',...
                slopeFitted,medianSlopeBoot,lowSlopeBoot,highSlopeBoot,round(100*options.bootConfInterval))});
        else
            title({sprintf('Threshold: %0.4f',thresholdFitted); sprintf('Slope: %0.4f',slopeFitted)})
        end
    end
    
    %% Get QuestPlus prediction and add to plot.
    %
    % Calculate QuestPlus prediction here and plot it.
    if ~isempty(options.questPara)
        predictedQuestPlus = qpPFWeibullLog(fineStimLevelsPlot',options.questPara);
        h_quest = plot(fineStimLevelsPlot,predictedQuestPlus(:,2),'k--','LineWidth',3);
    end
    
    % Add legend if you want.
    if (options.nBootstraps > 0)
        if ~isempty(options.questPara)
            legendHandles = [h_data h_pffit h_thresh h_bsthresh(2) h_bsthresh(1) h_quest];
        else
            legendHandles = [h_data h_pffit h_thresh h_bsthresh(2) h_bsthresh(1)];
        end
        if (options.addLegend)
            if ~isempty(options.questPara)
                legend(legendHandles, 'Data','PF-fit','PF-Threshold','BS-Threshold','BS-ConfInt', 'Quest-fit',...
                    'FontSize', 12, 'location', 'southeast');
            else
                legend(legendHandles, 'Data','PF-fit','PF-Threshold','BS-Threshold','BS-ConfInt', ...
                    'FontSize', 12, 'location', 'southeast');
            end
        end
    else
        if ~isempty(options.questPara)
            legendHandles = [h_data h_pffit h_thresh h_quest];
        else
            legendHandles = [h_data h_pffit h_thresh];
        end
        if (options.addLegend)
            if ~isempty(options.questPara)
                legend(legendHandles, 'Data','PF-fit','PF-Threshold','Quest-fit',...
                    'FontSize', 12, 'location', 'southeast');
            else
                legend(legendHandles, 'Data','PF-fit','PF-Threshold', ...
                    'FontSize', 12, 'location', 'southeast');
            end
        end
    end
    
    drawnow;
else
    % Print out the legendHandles empty when the 'verbose' is set to false.
    legendHandles = [];
end
