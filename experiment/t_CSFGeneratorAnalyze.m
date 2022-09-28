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

%% Start over.
clear; close all;

%% Set parameters here.
VERBOSE = true;
CHECKADAPTIVEMODE = false;
PF = 'weibull';

subjectName = 'Semin';
sineFreqCyclesPerDeg = 3;

%% Load the data and PF fitting.
%
% You can loop it if you want to fit multiple data. nData decides the
% number of data to load from the most recent one.
%
% Set startData to 0 if you want to read the data from the most recent.
startData = 0;
nData = 1;
SUBPLOT = false;
sizeSubplot = [round(nData/2) 4];

for dd = 1:nData
    % Load the data.
    if (ispref('SpatioSpectralStimulator','TestDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','TestDataFolder');
        testFilename = GetMostRecentFileName(fullfile(testFiledir,subjectName),...
            sprintf('CS_%s_%d_cpd',subjectName,sineFreqCyclesPerDeg),'olderDate',startData+dd-1);
        theData = load(testFilename);
    else
        error('Cannot find data file');
    end
    
    % Pull out the data here.
    nTrials = theData.estimator.nRepeat;
    [stimVec, responseVec, structVec] = combineData(theData.estimator);
    
    % Here we check if the Adaptive mode works fine if you want.
    if (CHECKADAPTIVEMODE)
        if (SUBPLOT)
            subplot(sizeSubplot(1),sizeSubplot(2),2*dd-1); hold on;
        end
        figure; clf;
        % Set variables here.
        nTrial = theData.estimator.nTrial;
        testTrials = linspace(1,nTrial,nTrial);
        testContrasts = 10.^stimVec;
        testPerformances = {structVec.outcome};
        
        % Set the color of the data point according to the subject's
        % performance. We will differenciate the marker point based on
        % each response either correct of incorrect.
        %
        % In testPerformances, correct response is allocated to 2 and
        % incorrect to 1. We will differenciate the marekr color based on
        % these.
        responseCorrect   = 2;
        responseIncorrect = 1;
        markerColorCorrect   = 5;
        markerColorIncorrect = 10;
        
        for tt = 1:nTrial           
            markerFaceColor(1,tt) = testPerformances{tt};
            markerFaceColor(find(markerFaceColor == responseCorrect)) = markerColorCorrect;
            markerFaceColor(find(markerFaceColor == responseIncorrect)) = markerColorIncorrect;
        end
        
        % Plot it.
        markerSize = 40;
        scatter(testTrials, testContrasts, markerSize, markerFaceColor, 'filled', 'MarkerEdgeColor', zeros(1,3));
        xlabel('Number of Trials', 'fontsize', 15);
        ylabel('Contrast', 'fontsize', 15);
        title('Checking Adaptive Mode', 'fontsize', 15);
        legend('Blue = correct / Yellow = incorrect','location','northeast','fontsize',15);
        
        clear markerFaceColor;
    end
    
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
    
    thresholdCriterion = 0.81606;
    [threshold, para, dataOut] = theData.estimator.thresholdMLE(...
        'thresholdCriterion', thresholdCriterion, 'returnData', true);
    
    % Set the contrast levels in linear unit.
    examinedContrastsLinear = 10.^dataOut.examinedContrasts;
    
    % PF fitting here.
    if (SUBPLOT)
        subplot(sizeSubplot(1),sizeSubplot(2),dd*2); hold on;
    end
    [paramsFitted(:,dd)] = FitPFToData(examinedContrastsLinear, dataOut.pCorrect, ...
        'PF', PF, 'nTrials', nTrials, 'verbose', VERBOSE, 'figureWindow', ~SUBPLOT, 'pointSize', pointSize);
    
    clear pointSize;
end

% Save the graph if you want.
SAVEFITTING = false;

if (SAVEFITTING)
    fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Experiment'; 
    cd(fileDir);
    fileName = sprintf('%d_deg.tiff', sineFreqCyclesPerDeg);
    saveas(gcf, fileName);
end

%% Plot the CSF curve here.
CSFCURVE = false;

if (CSFCURVE)
    % Set target spatial frequency and threshold values. We type manually
    % here for now, but we can elaborate it later.
    
    % Semin
    testSpatialFrequency = [3 6 9 12 18];
    logTestSpatialFrequency = log10(testSpatialFrequency);
%     threshold_Semin = [0.0017 0.0032 0.0040 0.0043 0.0046 0.0072];
    threshold_Semin = [0.0031 0.0055 0.0059 0.0053 0.0139];
    sensitivity_Semin = log10(1./threshold_Semin);
    
%     threshold_David = [0.0015 0.0033 0.0043 0.0065];
    threshold_David = [0.0027 0.0050 0.0062 0.0052 0];
    sensitivity_David = log10(1./threshold_David);
    
    % Plot it.
    figure; clf; hold on;
    plot(logTestSpatialFrequency, sensitivity_Semin, 'g.-','markersize',20,'linewidth',2);
    plot(logTestSpatialFrequency, sensitivity_David, 'b.-','markersize',20,'linewidth',2);
    xlabel('Spatial Frequency (cpd)','fontsize',15);
    ylabel('Contrast Sensitivity','fontsize',15); 
    xticks(logTestSpatialFrequency);
    xticklabels(testSpatialFrequency);
    yticks(sort([sensitivity_David sensitivity_Semin]));
    yticklabels(round(10.^sort([sensitivity_David sensitivity_Semin])));
    legend('Semin','David','fontsize',15);
end

% Save the graph if you want.
SAVEFITTING = false;
if (SAVEFITTING)
    fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Experiment'; 
    cd(fileDir);
    fileName = sprintf('CSF_curves.tiff', sineFreqCyclesPerDeg);
    saveas(gcf, fileName);
end