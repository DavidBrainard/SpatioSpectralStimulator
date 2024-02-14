function [f0] = FindInitialFrequencyToFitSineWave(waveform,options)
% This function searches initial frequency to fit the sine waveform.
%
% Syntax:
%    [f0] = FindInitialFrequencyToFitSineWave(waveform)
%
% Description:
%    This function runs a loop to grid search an initial frequency value to
%    fit sinusoidal waveform.
%
% Inputs:
%     waveform                  - Target signal you want to fit.
%
% Outputs:
%     f0                        - Found initial frequency value to fit the
%                                 sine to the given waveform.
%
% Optional key/value pairs:
%    SF                         - Spatial frequency of the waveform to fit.
%                                 This will decide the initial set f0 value
%                                 of the searching loop. Default to 3 (cpd).
%    f0_searchInterval          - Interval size within the grid search of
%                                 the loop. Default to 0.05.
%    verbose                    - Display additional information (default: false).

% History:
%    12/11/23  smo             - Wrote it.

%% Set variables.
arguments
    waveform
    options.SF (1,1) = 3
    options.f0_searchInterval (1,1) = 0.0025
    options.verbose (1,1) = true
end

%% Set the waveform class as double.
if ~strcmp(class(waveform),'double')
    waveform = double(waveform);
end

%% Here we runs a loop to grid search an initial frequency value.
%
% Set the initial f0 value differently over spatial frequency. We can
% always start from 0, but setting this value differently (higher value)
% will speed up the searching process.
switch options.SF
    case 1 
        f0_lb = 0;
    case 3
        f0_lb = 0;
    case 6
        f0_lb = 3;
    case 9
        f0_lb = 4;
    case 12
        f0_lb = 8;
    case 18
        f0_lb = 20;
    otherwise
        f0_lb = 0;
end

% Make a loop to search from here. We will use f0 value one by one to fit
% the sine wave until it finds the curve that has the certain level of
% correaltion.
f0 = f0_lb;
trial = 1;
while 1
    % Fitting happens here.
    [params, fittedSignal] = FitSineWave(waveform,'f0',f0,'verbose',false,'FFT',false);
    
    % Set the target correlation differently. For lower spatial frequency,
    % it's not possible to achieve as high as 99% correlation between two
    % curves.
    switch options.SF
        case 1
            targetCorrSignals = 0.92;
        case 3
            targetCorrSignals = 0.94;
            % For SACCSFA185.
%             targetCorrSignals = 0.941; 
        case 6
            targetCorrSignals = 0.94;
%             targetCorrSignals = 0.937;
        case 9
            targetCorrSignals = 0.965;
            
            % For Camera MTF.
%             targetCorrSignals = 0.945;
        case 12
            targetCorrSignals = 0.98;
        case 18
            targetCorrSignals = 0.97;
        otherwise
            targetCorrSignals = 0.97;
    end
    
    % Stop the loop if we found a fit exceeding the criteria. We find
    % it based on the correlation between the original and fitted
    % signal.
    fittedCorrSignals = corr(waveform',fittedSignal');
    if fittedCorrSignals > targetCorrSignals
        break;
    end
    
    % If not, we will update f0 value and keep searching.
    f0 = f0 + options.f0_searchInterval;
    trial = trial+1;
    
    % Show progress every 10 trials.
    if (options.verbose)
        if mod(trial,10) == 0
            fprintf('Fitting progress - Number of trials = (%d) \n',trial);
        end
    end
end

% Plot the results.
if (options.verbose)
    figure; hold on;
    plot(waveform,'b-');
    plot(fittedSignal,'r-');
    title(sprintf('Found f0 is (%.4f)',f0));
    legend('Origianl','Fit');
    ylim([0 max(waveform)*1.05]);
end

end
