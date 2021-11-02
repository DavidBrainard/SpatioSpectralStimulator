function [spdMeasured] = MeasureProjectorPlainScreenSettings(theSettings,S,window,windowRect,options)
% Measure the SPD over projector primary settings.
%
% Syntax: [[spdMeasured] = MeasureProjectorPlainScreenSettings(theSettings,S,window,windowRect)
%
% Description:
%    This measures the SPD according to projector primary settings. In case
%    you use multi-subprimary device, you should set subprimary settings in
%    advance before calling this function.
%
% Inputs:
%    theSettings -                Projector input settings to measure for.
%    S -                          Spectrum measurement range in wavelength.
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Outputs:
%    spdMeasured -                Measurement results of the SPDs for the
%                                 contrast testing points.
%
% Optional key/value pairs:
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.

% History:
%    10/08/21  smo                Started on it
%    10/14/21  smo                Made a working draft.
%    10/15/21  smo                Fixed the error, clean it, better naming
%                                 the variables.
%    10/18/21  smo                Made it simpler and it takes projector
%                                 settings and save out the spd data only.
%    10/19/21  smo                Now it does not take subprimary settings
%                                 as input. It only takes care of the
%                                 projector primary settings.
%    10/26/21  smo                Delete calstruct data in input, and add
%                                 'S' as input for getting spectrum
%                                 measurement range.
%    11/02/21  smo                Clean it after deleting old features.

%% Set parameters.
arguments
    theSettings
    S
    options.measurementOption (1,1) = true
    options.verbose (1,1) = true
end

%% Set the projector mode.
%
if (options.measurementOption)
    if (options.projectorMode)
        commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Normal mode (Default)
        unix(commandNormal)
        disp('Projector is set as Normal mode');
    else
        commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit'; % Steady-on mode
        unix(commandSteadyOn)
        disp('Projector is set as Steady-on mode');
    end
end

% Number of contrast test points.
nTestPoints = size(theSettings,2);

%% Display the projector image and measure it.
% Set projector setting for each test point on the gabor patch and measure it.
if (options.measurementOption)
    % Measure the test points.
    for tt = 1:nTestPoints
        % Set the projector settings and display it as a plane screen.
        SetProjectorPlainScreenSettings(theSettings(:,tt),window,windowRect,'verbose',options.verbose); 

        % Measure it.
        spdMeasured(:,tt) = MeasureSPD('S',S);
        if (options.verbose)
            fprintf('           Measurement complete! - Test Point (%d/%d) \n',tt,nTestPoints);
        end
    end
else
    % Just record zero spectra for both test and background
    % when skipping the measurements.
    spdMeasured = zeros(S(3),nTestPoints);
    if (options.verbose)
        fprintf('           Measurement has been skipped \n');
    end
end

%% Plot the results.
if (options.verbose)
    figure; hold on
    wls = SToWls(S);                          % Spectrum range.
    plot(wls,spdMeasured,'r-','LineWidth',2); % Test points spds.
    title('Measured SPDs');
    xlabel('Wavelength (nm)')
    ylabel('Spectral Intensity');
end
end