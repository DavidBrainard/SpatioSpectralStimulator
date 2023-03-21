function [fitError] = SmoothnessSearchErrorFunction(smoothingParam, csfSFs, fitData, evaluateData)
% It computes the cross validated error when fitting smoothing spline
% function.
%
% Syntax:
%    [fitError] = SmoothnessSearchErrorFunction(x, csfSFs, fitData, evaluateData)
%
% Description:
%    This is computing the fitting error in the same way that we did in our
%    grid search over smoothing parameter values. It is in the form of
%    function so that fmincon understands.
%
% Inputs:
%    smoothingParam
%    csfSFs
%    fitData
%    evaluateData
%
% Outputs:
%    fitError
%
% Optional key/value pairs.
%    verbose
%
% See also:
%    SACC_FitCSF.

% History:
%   03/17/23  dhb, smo         - Started on it.
%   03/20/23  smo              - Made it as a separate function.

%% Get the number of spatial frequency and data to fit.
nSFs = length(csfSFs);
nFits = size(fitData,2);

%% Fitting happens here.
%
% Use the spline to fit the fit data, for each cross valiation iteration.
% Each fit/evaluate csf pair in a column of fitData/evaluteData.
fitError = 0;
for ii = 1:nFits
    % Fit here and get predicted data.
    smoothFit = fit(csfSFs',fitData{ii}','smoothingspline','SmoothingParam',smoothingParam);
    prediction{ii} = feval(smoothFit,csfSFs);
    
    % Evaluate against matched evaluation data.
    diff = evaluateData{ii}' - prediction{ii};
    fitError = fitError + sum(diff.^2);
end

%% Normalize the fitting error.
fitError = sqrt(fitError/(nFits * nSFs));

end