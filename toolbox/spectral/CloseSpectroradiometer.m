function [] = CloseSpectroradiometer(options)
% This closes the spectroradiometer after finishing the measurements.
% Syntax: [[spdMeasured] = MeasurePlainScreenSettings(theSettings,S,window,windowRect)
%
% Description:
%    This simply disconnects your spectroradiometer. This should be used
%    after finishing up your measurements so that you don't need to connect
%    the device continuosuly when you repeat your measurements.
%
% Inputs:
%    N/A
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    'portNum' -                  Set port number that your device was
%                                 connected to. Default to 5 (PR670).
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% See also:
%    N/A

% History:
%    11/05/21  smo                Started on it

%% Set parameters
arguments
   options.portNum (1,1) = 5 
   options.verbose (1,1) = true
end

%% Set parameters.
CMClose(options.portNum);

if (options.verbose)
    switch (options.portNum)
        case 5
             disp('Successfully disconnected to PR-670!')
    end
end

end