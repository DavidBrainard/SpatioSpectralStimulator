function [measuredSpd] = MeasureSpectroradiometer(options)
% Measure the SPD using the spectroradiometer.
%
% Syntax:
%    [measuredSpd] = MeasureSpectroradiometer()
%
% Description:
%    This measures the spectral power distribution using the
%    spectroradiometer. This function should be used whenever measurement
%    is needed in SACC project.
%
% Inputs:
%    N/A
%
% Outputs:
%    measuredSpd -                Measured spd result for the given
%                                 wavelength range.
%
% Optional key/value pairs:
%    'S' -                        Measurement wavelength range. Default
%                                 [380 2 201], which means measuring from
%                                 380 nm with 2 nm interval to have a total
%                                 of 201 components.
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%    'portNum' -                  Port number connected to the
%                                 spectroradioment using the function
%                                 CMCeckInit. Default to 5 (PR670).
%    'verbose'           -        Boolean. Default true.  Controls plotting
%                                 and printout.
%
% See also: 
%    OpenRadiometer.

% History:
%    11/01/21  smo                Started on it.
%    11/02/21  smo                Separate the pre-measure part and
%                                 measurement as separate functions.

%% Set parameters.
arguments
    options.S (1,3) = [380 2 201]
    options.measurementOption (1,1) = true
    options.portNum (1,1) = 5
    options.verbose (1,1) = true
end

%% Measure it.
if (options.measurementOption)
    measuredSpd = MeasSpd(options.S, options.portNum, 'all');
  
if (options.verbose)
        disp('Measurement complete!');
    end
else
    measuredSpd = zeros(options.S(3),1);
    if (options.verbose)
        disp('Measurement has been skipped!');
    end
end
end
