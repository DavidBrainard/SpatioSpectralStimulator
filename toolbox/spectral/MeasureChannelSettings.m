function [spdMeasured] = MeasureChannelSettings(theSettings,S,options)
% Measure the spectrum for passed channel settings.
%
% Syntax:
%    [spdMeasured] = MeasureChannelSettings(theSettings,S)
%
% Description:
%    This measures the spd that we actually get.  Assumes that the plain
%    screen has been set up in advance, as well as the projector mode set.
%
% Inputs:
%    theSettings -                Channel settings we wish to set.  This
%                                 is a matrix with nPrimaries columns (one for
%                                 each screen primary) and nChannels
%                                 rows, one for each channel primary. These
%                                 values are specified as real numbers
%                                 between 0 and 1, and are assumed to be
%                                 after gamma correction.
%    S -                          Spectrum measurement range in wavelength.
%
% Outputs:
%    targetSpdMeasured -          Measured SPDs of the desired isolating
%                                 target primaries.
%
% Optional key/value pairs:
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%    'nPrimaries -                The number of primaries in the displaying
%                                 device.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% See also:
%    SetChannelSettings.

% History:
%    11/16/21  dhb, smo     Wrote it.

%% Set parameters.
arguments
    theSettings {mustBeInRange(targetPrimaries,0,1,"inclusive")}
    S (1,3)
    options.measurementOption (1,1) = true
    options.nPrimaries (1,1) = 3
    options.verbose (1,1) = true
end

%% Set primary settings and measure it.
if (options.measurementOption)
    % Set primary input settings. This measures one primary at a time, so
    % the other primaries are set to zero.
    SetChannelSettings(theSettings);
    spdMeasured = MeasureSpectroradiometer('S',S);
    if (options.verbose)
        fprintf('Made the measurement!\n');
    end
else
    % Save out the zero spectrum if the measurement is skipped.
    spdMeasured = zeros(S(3),1); 
    if (options.verbose)
        fprintf('Measurement has been skipped!\n');
    end
end

% Plot the results.
if (options.verbose)
    figure; hold on;
    plot(SToWls(S),spdMeasured,'r-','LineWidth',2); 
    xlabel('Wavelength (nm)');
    ylabel('Spectral power');
    legend('Target','Measurement');
    title(sprintf('Primary %d',targetPrimaryNum));
end

end