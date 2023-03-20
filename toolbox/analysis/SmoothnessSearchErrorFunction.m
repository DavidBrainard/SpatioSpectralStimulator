% This function computes the cross validated error.  This is really just
% computing the same thing you are already plotting in your grid search
% over smoothness parameter values, placed into a function that fmincon
% understands.
function [fitError] = SmoothnessSearchErrorFunction(x, csfSFs, fitData, evaluateData)
% Grab current smoothness param from parameter vector.
smoothingParam = x(1);

% Use the spline to fit the fit data, for each cross
% valiation iteration Each fit/evaluate csf pair in a
% column of fitData/evaluteData.
nSFs = length(csfSFs);
nFits = size(fitData,2);
fitError = 0;

% Use spline to fit the data with passed smoothing param.
for ii = 1:nFits
    smoothFit = fit(csfSFs',fitData{ii}','smoothingspline','SmoothingParam',smoothingParam);
    prediction{ii} = feval(smoothFit,csfSFs);
    
    % Evaluate against matched evaluation data.
    diff = evaluateData{ii}' - prediction{ii};
    fitError = fitError + sum(diff.^2);
end

% Normalize the error.
fitError = sqrt(fitError/(nFits * nSFs));

end