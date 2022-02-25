function [paramsFitted] = FitPFToData(stimLevels,pCorrect,options)
% Fit Psychometric function to the given data.
%
% Syntax:
%    [] = FitPFToData()
%
% Description:
%    This fits Psychometric function to the data given as an input.
%
% Inputs:
%    stimLevels -                 ddd
%    pCorrect - 
%
% Outputs:
%    paramsFitted
%
% Optional key/value pairs:
%    PF -                         Default to 'weibull'. Choose the function
%                                 to fit the data either 'weibull' or
%                                 'logistic'.
%    verbose -                    Boolean. Default true.  Controls plotting
%                                 and printout.
%    nTrials -                    Default to 20. The number of trials
%                                 conducted per each stimulus level.
%
% See also:
%    N/A

% History:
%   02/25/22 dhb, smo             Started on it.

%% Set parameters.
arguments
    stimLevels
    pCorrect
    options.PF = 'weibull'
    options.nTrials (1,1) = 20
    options.thresholdCriterion (1,1) = 0.81606
    options.verbose (1,1) = true
end

%% Set fitting type.
%
% Choose a function to use.
switch options.PF
    case 'weibull'
        PF = @PAL_Weibull;
    case 'logistic'
        PF = @PAL_Logistic
    otherwise
end

%% Fitting PF here.
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

% Set which params to be found.
% paramsFree = [thresh slope guess lapse]; 0 = fixed; 1 = free.
paramsFree = [1 1 0 1]; 
paramsFitted = PAL_PFML_Fit(stimLevels, nCorrect, ...
    nTrialsPerContrast, searchGrid, paramsFree, PF);

% Make a smooth curves with finer stimulus levels.
nFineStimLevels = 1000;
fineStimLevels = linspace(0,max(stimLevels),nFineStimLevels);
smoothPsychometric = PF(paramsFitted, fineStimLevels);
thresholdFitted = PF(paramsFitted, options.thresholdCriterion, 'inv');

%% Plot the results.
%
% All data.
figure; clf; hold on;
marekrColorGray = [0.7 0.7 0.7];
plot(stimLevels,pCorrect,'ko','MarkerFaceColor',marekrColorGray,'MarkerSize',10);
plot(fineStimLevels,smoothPsychometric,'r','LineWidth',3);

% Threshold point.
plot(thresholdFitted,options.thresholdCriterion,'ko','MarkerFaceColor','r','MarkerSize',12);
ylim([0 1]);
xlabel('Contrast');
ylabel('pCorrect');
legend('Data','PF fit','Threshold');

end
