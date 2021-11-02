function [measuredSpd] = MeasureSPD(options)
% Measure the SPD using the spectroradiometer.
%
% Syntax:
%    [measuredSpd] = MeasureSPD()
%
% Description:
%    This measures the spectral power distribution using the
%    spectroradiometer. This code connects to PR670, and this function
%    should be used whenever measurement is needed in SACC project.
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
%    'verbose'           -        Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    11/1/21  smo                 Started on it.

%% Set parameters.
arguments
    options.S (1,3) = [380 2 201]
    options.measurementOption (1,1) = true
    options.verbose (1,1) = true
end

%% Connect to the spectroradiometer PR670.
% Give read-and-write permission to the spectroradiometer.
% This part is exclusivley when using the Linux box, and
% the PR670 is connected to the device either 'ttyACM0' or 'ttyACM1',
% so give permission to both here.
if (options.measurementOption)
    commandConnectToPR670_1 = 'sudo chmod a+rw /dev/ttyACM0';
    unix(commandConnectToPR670_1)
    commandConnectToPR670_2 = 'sudo chmod a+rw /dev/ttyACM1';
    unix(commandConnectToPR670_2)
    
    % Connect to spectroradiometer PR670.
    portNumPR670 = 5;
    CMCheckInit(portNumPR670);
    
    % Measure it.
    measuredSpd = MeasSpd(options.S,portNumPR670,'all');
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
