function [params fittedWaveform] = FitSineWave(waveform,options)
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
%     waveform                  - Target signal you want to fit.
%
% Outputs:
%    params                     - Parameters fitted. It includes three
%                                 params, A (amplitute), f (frequency), phi
%                                 (phase shift).
%
% Optional key/value pairs:
%    A0                        - Initial amplitude (default: 1).
%    f0                        - Initial frequency (default: 1).
%    phi0                      - Initial phase (default: 0).
%    lb                        - Lower bounds for parameters [A, f, phi] (default: [0, 0, -pi]).
%    ub                        - Upper bounds for parameters [A, f, phi] (default: [Inf, Inf, pi]).
%    verbose                   - Display additional information (default: false).

% History:
%    11/09/23  smo              - Started on it

%% Set variables.
arguments
    waveform
    options.A0 (1,1) = 1
    options.f0 (1,1) = 1
    options.phi0 (1,1) = 0
    options.lb (1,3) = [0, 0, -pi];
    options.ub (1,3) = [Inf, Inf, pi];
    options.verbose (1,1) = false
end

%% Fourier transform.
samplingRate = length(waveform);
duration = 1;
t = linspace(0,duration,samplingRate);

% Perform the Fourier Transform here.
N = length(waveform);
frequencies = (0:N-1) * (samplingRate / N);
signal_fft = fft(waveform);
signal_fft_magnitude = abs(signal_fft) / N;

% Find the fundamental frequencies
[~, idxFundamental] = max(signal_fft_magnitude);
fundamental_frequency = frequencies(idxFundamental);

% Plot the original signal and its frequency spectrum
if (options.verbose)
    figure;
    xlim([-1 N]);
    subplot(2, 1, 1);
    plot(t, signal);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Target signal');
    
    subplot(2, 1, 2);
    plot(frequencies, signal_fft_magnitude);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Frequency Spectrum');
    
    fprintf('Fundamental Frequency 1: %.2f Hz\n', fundamental_frequency);
end

% Get the reconstructrued signal.
reconstructed_signal1 = sin(2 * pi * fundamental_frequency * t);

% Plot the original and reconstructed signals
figure;
subplot(2, 1, 1);
plot(t, waveform);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');

subplot(2, 1, 2);
plot(t, reconstructed_signal1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal (Fundamental Frequency 1)');

% Fit the sinusoidal signal using fmincon from here.
%
% Define the observed waveform data (replace this with your actual data)
x_data = t;
observedWaveform = waveform;

% Define the objective function for optimization
objective_function = @(params) norm(params(1)*sin(2*pi*params(2)*x_data + params(3)) - observedWaveform);

% Set the initial starting value.
initial_guess = [options.A0, options.f0, options.phi0];

% Call fmincon to optimize the parameters
optionsFmincon = optimoptions('fmincon', 'Display', 'iter');
optimized_params = fmincon(objective_function, initial_guess, [], [], [], [], options.lb, options.ub, [], optionsFmincon);

% Extract optimized parameters
A_optimized = optimized_params(1);
f_optimized = optimized_params(2);
phi_optimized = optimized_params(3);

% Save out the fitted params.
params(1) = A_optimized;
params(2) = f_optimized;
params(3) = phi_optimized;

% Generate the fitted waveform using the optimized parameters
fittedWaveform = A_optimized * sin(2*pi*f_optimized*x_data + phi_optimized);

if (options.verbose)
    % Plot the observed and fitted waveforms
    figure;
    plot(x_data, observedWaveform, 'b', 'DisplayName', 'Observed Waveform');
    hold on;
    plot(x_data, fittedWaveform, 'r', 'DisplayName', 'Fitted Waveform');
    xlabel('x');
    ylabel('Amplitude');
    legend('show');
    grid on;
    
    % Display optimized parameters
    fprintf('Optimized Parameters:\n');
    fprintf('a = %.4f\n', A_optimized);
    fprintf('b = %.4f\n', f_optimized);
    fprintf('c = %.4f\n', phi_optimized);
end

end