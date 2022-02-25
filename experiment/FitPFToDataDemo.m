    clear; close all;

    PF = @PAL_Weibull;
    simParams = [0.02 4 0.5 0.01];

    contrasts = linspace(0.01,0.05,10);
    pCorrect = PF(simParams,contrasts);
    
    figure; clf; hold on
    plot(contrasts,pCorrect,'ro');


    % Set up the PF fitting (requires Palamedes toolbox).  Note that the
    % catch trials are added in here in the call to the fit.
    nTrials = 20;
    nTrialsPerContrast = nTrials*ones(size(contrasts));
    nCorrect = round(pCorrect.* nTrialsPerContrast);

    
    searchGrid = [mean(contrasts) 2 0.5 0.01];

    % Fit
    paramsFree = [1 1 0 0]; %[thresh slope guess lapse]; 0 = fixed; 1 = free
        paramsFitted = PAL_PFML_Fit(contrasts, nCorrect, ...
            nTrialsPerContrast, searchGrid, paramsFree, PF);

    nFineContrasts = 1000;
    fineContrasts = linspace(min(contrasts),max(contrasts),nFineContrasts);
    smoothPsychometric = PF(paramsFitted, fineContrasts);

    thresholdCriterion = 0.82;
    threshold = PF(paramsFitted, thresholdCriterion, 'inv');

    plot(fineContrasts,smoothPsychometric,'k');
    plot(threshold,thresholdCriterion,'ko','MarkerFaceColor','k','MarkerSize',12);

   