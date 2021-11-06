function [] = OpenSpectroradiometer(options)
% This makes the radiometer ready to use.
%
% Syntax:
%    [] = OpenSpectroradiometer()
%
% Description:
%    This makes the radiometer ready to use. It connects to PR670
%    (default), and this should be used whenever measurement is
%    needed in SACC project.
%
% Inputs:
%    N/A
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%    'portNum' -                  Port number connected to the
%                                 spectroradioment using the function
%                                 CMCeckInit. Default to 5 (PR670).
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% See also: 
%    MeasureSpd.

% History:
%    11/2/21  smo                 Started on it.

%% Set parameters.
arguments
    options.measurementOption (1,1) = true
    options.portNum (1,1) = 5
    options.verbose (1,1) = true
end

%% Connect to the spectroradiometer.
% Give read-and-write permission to the spectroradiometer. This part is
% exclusivley when using the Linux box, and the PR670 is connected to the
% device either 'ttyACM0' or 'ttyACM1', so give permission to both here.
% And you need to type the account password (colorlab) in the command space
% to do so.
if (options.measurementOption)
    commandConnectToPR670_dev1 = 'sudo chmod a+rw /dev/ttyACM0';
    connection_dev1 = unix(commandConnectToPR670_dev1);
    commandConnectToPR670_dev2 = 'sudo chmod a+rw /dev/ttyACM1';
    connection_dev2 = unix(commandConnectToPR670_dev2);
    
    % Check which dev to connect to PR670.
    if connection_dev1 == 1
       if (options.verbose)
          disp('Spectroradiometer is connected to ttyACM0'); 
       end
    elseif connection_dev2 == 1
       if (options.verbose)
          disp('Spectroradiometer is connected to ttyACM1');
       end
    else
        error('Spectroradiometer is not found!');
    end
    
    % Connect to spectroradiometer. Default to 5 (PR670).
    CMCheckInit(options.portNum);
end
end