function [paramsFitted] = FitPFToData(stimLevels,pCorrect,options)
% Fit Psychometric function to the given data.
%
% Syntax:
%    [paramsFitted] = FitPFToData(stimLevels,pCorrect)
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
%    nTrials -                    Default to 20. The number of trials
%                                 conducted per each stimulus level.
%    thresholdCriterion -         Default to 0.81606. This is the value of
%                                 pCorrect as a criteria to find threshold.
%    figureWindow -               Default to true. Make a new figure window
%                                 to plot the results. If you want to plot
%                                 multiple threshold results in a subplot,
%                                 set this to false and do so.
%    verbose -                    Default to true. Boolean. Controls
%                                 plotting and printout.
%
% See also:
%    N/A

% History:
%   02/25/22 dhb, smo             Started on it.
%   03/14/22 smo                  Added a plotting option to make different
%                                 marker size over different number of
%                                 trials.

%% Set parameters.
arguments
    stimLevels
    pCorrect
    options.PF = 'weibull'
    options.paramsFree (1,4) = [1 1 0 1]
    options.nTrials (1,1) = 20
    options.thresholdCriterion (1,1) = 0.81606
    options.figureWindow (1,1) = true
    options.verbose (1,1) = true
    options.pointSize = ones(1,length(stimLevels))*100
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
        PF = @PAL_Logistic
    otherwise
end

%% Fitting PF.
%
% Set up the PF fitting (requires Palamedes toolbox).  Note that the
% catch trials are added in here in the call to the fit.
nTrialsPerContrast = options.nTrials * ones(size(stimLevels));
nCorrect = round(pCorrect .* nTrialsPerContrast);

% Set an initinal search parameters.
threshold = mean(stimLevels);
slope = 2;
guess = 0.5;
lapse = 0.01;
searchGrid = [threshold slope guess lapse];

% PF fitting happens here.
paramsFitted = PAL_PFML_Fit(stimLevels, nCorrect, ...
    nTrialsPerContrast, searchGrid, options.paramsFree, PF);

% Make a smooth curves with finer stimulus levels.
nFineStimLevels = 1000;
fineStimLevels = linspace(0, max(stimLevels), nFineStimLevels);
smoothPsychometric = PF(paramsFitted, fineStimLevels);
thresholdFitted = PF(paramsFitted, options.thresholdCriterion, 'inv');

if (options.verbose)
    fprintf('Threshold was found at %.4f (linear unit) \n', thresholdFitted);
end

%% Plot the results if you want.
if (options.verbose)
    if (options.figureWindow)
        figure; clf; hold on;
    end 
    
    % Plot all experimental data (gray points).
    %
    % Marker size will be different over the number of the trials per each
    % test point.
    scatter(stimLevels, pCorrect, options.pointSize,...
        'MarkerEdgeColor', zeros(1,3), 'MarkerFaceColor', ones(1,3) * 0.5, 'MarkerFaceAlpha', 0.5);
    plot(fineStimLevels,smoothPsychometric,'r','LineWidth',3);
    
    % Mark the threshold point (red point).
    plot(thresholdFitted,options.thresholdCriterion,'ko','MarkerFaceColor','r','MarkerSize',12);
    ylim([0 1]);
    xlabel('Contrast', 'FontSize', 15);
    ylabel('pCorrect', 'FontSize', 15);
    legend('Data','PF fit','Threshold','FontSize', 12, 'location', 'southeast');
    title(append('Threshold: ', num2str(round(thresholdFitted,4))), 'FontSize', 15);
end

end
