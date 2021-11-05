function [targetSpdMeasured] = MeasureDesiredTargetPrimaries(targetPrimaries,subPrimaryCalStructData, ...
    targetPrimaryNum,options)
% Measure the desired target primaries to use them for computing the contrast image.
%
% Syntax:
%    [targetSpdMeasuredRaw,targetSpdMeasuredNorm,measuredToTarget] = MeasureDesiredTargetPrimaries(targetPrimaries,subPrimaryCalStructData, ...
%                                                                    targetPrimaryNum,options)
%
% Description:
%    This measures the target primaries that we actually get and then use
%    the measured rather than the nominal primaries to compute the image.
%
% Inputs:
%    targetPrimaries -            Target primaries you wish to reproduce.
%                                 These have not yet been gamma corrected.
%    subPrimaryCalStructData -    The calibration object that describes the
%                                 device we're working with.
%    targetPrimaryNum -           Target primary number. This is required
%                                 when setting up the projector current input
%                                 levels.
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
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either to be 'Normal' (true) or
%                                 'Steady-on' (false).
%    'nPrimaries -                The number of primaries in the displaying
%                                 device.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% See also:
%    SetSubprimarySettings.

% History:
%    10/06/21  smo                Started on it
%    10/07/21  smo                Makes it working for a single primary and
%                                 added a feature to skip the measurement part.
%    10/19/21  smo                Added to save out the normalized spd results.
%    10/20/21  smo                Discard the normalized spd outputs.
%    11/02/21  smo                Clean it by deleting old features.

%% Set parameters.
arguments
    targetPrimaries {mustBeInRange(targetPrimaries,0,1,"inclusive")}
    subPrimaryCalStructData
    targetPrimaryNum (1,1)
    options.measurementOption (1,1) = true
    options.projectorMode (1,1) = true
    options.nPrimaries (1,1) = 3
    options.verbose (1,1) = true
end

%% Get number of discrete input levels for the device.
% So, control values go from 0 to (subprimaryNInputLevels-1).
subprimaryNInputLevels = size(subPrimaryCalStructData.get('gammaInput'),1);

% Measurement range.
S = subPrimaryCalStructData.get('S');

% Set primary and subprimary numbers.
nSubprimaries = subPrimaryCalStructData.get('nDevices');
if (targetPrimaryNum > 3) || (targetPrimaryNum < 1)
    error('Target primary number should be set within [1, 2, 3]')
end

% Target Spd.
% This will be compared to the measurement result at the end.
targetSpd = PrimaryToSpd(subPrimaryCalStructData,targetPrimaries);

%% Set primary settings and measure it.
if (options.measurementOption)
    % Set primary input settings. This measures one primary at a time, so
    % the other primaries are set to zero.
    targetSubprimarySettings = zeros(nSubprimaries,options.nPrimaries);
    targetSubprimarySettings(:,targetPrimaryNum) = PrimaryToSettings(subPrimaryCalStructData,targetPrimaries);
    SetSubprimarySettings(targetSubprimarySettings,'projectorMode',options.projectorMode);
    
    % Measure it.
    if (options.verbose)
        for ss = 1:nSubprimaries
            fprintf('Subprimary settings %2.0f - %.2f\n',ss,targetSubprimarySettings(ss));
        end
    end
    targetSpdMeasured = MeasureSpectroradiometer('S',S);
    if (options.verbose)
        fprintf('Measurement complete! - Primary %d\n',targetPrimaryNum);
    end
else
    % Save out the zero spectrum if the measurement is skipped.
    targetSpdMeasured = zeros(S(3),1); 
    if (options.verbose)
        fprintf('Measurement has been skipped!\n');
    end
end

% Plot the results.
if (options.verbose)
    figure; hold on;
    plot(SToWls(S),targetSpd,'k-','LineWidth',3); % Target Spd.
    plot(SToWls(S),targetSpdMeasured,'r-','LineWidth',2); % Measured raw spd.
    xlabel('Wavelength (nm)');
    ylabel('Spectral power');
    legend('Target','Measurement_Raw','Meausurement_Norm');
    title(sprintf('Primary %d',targetPrimaryNum));
end
end