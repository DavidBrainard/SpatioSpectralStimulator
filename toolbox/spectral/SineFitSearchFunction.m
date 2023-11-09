function [fitError] = SineFitSearchFunction(signal,t,x)
% A form of function to fit sine signal.
%
% Syntax:
%    [params] = FitSineSignal(signal)
%
% Description:
%    This is to fit sinusoidal waveform. We wrote this to characterize the
%    chromatic aberration of the SACCSFA system for both LCA and TCA.
%
% Inputs:
%     signal                    - Target signal you want to fit.
%
% Outputs:
%    params                     - Parameters fitted. It includes three
%                                 params, A (amplitute), f (frequency), phi
%                                 (spatial shift).
%
% Optional key/value pairs:
%

% History:
%    11/09/23  smo              - Started on it

%% Set variables.
% arguments
%     signal
%     x (1,3)
%     fitError
% end

%% Set the form of the fitting function.
%
% We will set the initial value for each parameter.
fittedWaveform = x(1) * sin(2 * pi * x(2) * t + x(3));

%% Calculate the fit error.
fitError = sum(signal-fittedWaveform);

end
