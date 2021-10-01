function [] = MeasureDesiredTargetPrimaries(foo,projectorMode)
% Measure the desired target primaries to use them for computing the image.
% 
% Syntax:
%    foo
%
% Description:
%    foo
%
% Inputs:
%    foo -                        To be filled
%
% Outputs:
%    foo -                        To be filled
%
% Optional key/value pairs:
%    N/A
%
% History: 
%    09/29/21  smo  Started on it

%% This is where we would measure the primaries we actually get and then use
% the measured rather than the nominal primaries to compute the image.

%% Initialize (connect to both display and measurement device).
%
% This part prepares for the measurements which includes connecting PC to the
% projector, setting the projector display state (normal/steady-on modes),
% and talking to the measurement device to ready to run (PR670).
%
% Add VPixx toolbox 'Datapixx' to the path.
toolboxDirectory = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx'; % This is where Linux box store the file.
addpath(genpath(toolboxDirectory)); 

% Connect to the Vpixx projector
isReady = Datapixx('open');
isReady = Datapixx('IsReady'); 

% Set the projector mode.
switch lower(projectorMode);
    case 'normal'
        commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Normal mode (Default)
        unix(commandNormal)
        disp('Projector is set as Normal mode');
    case 'steady'
        commandSteadyOn = 'vputil rw 0x1c8 0x7 -q quit'; % Steady-on mode
        unix(commandSteadyOn)
        disp('Projector is set as Steady-on mode');
    otherwise
        commandNormal = 'vputil rw 0x1c8 0x0 -q quit'; % Default setting as Normal mode
        unix(commandNormal)
        disp('Projector is set as Normal mode as default');
end

% Give read-and-write permission to the spectroradiometer.
% This part is exclusivley when using the Linux box, and 
% the PR670 is connected to either 'ttyACM0' or 'ttyACM1', so give
% permission to both here.
command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)
command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)
command_check = 'ls -l /dev/ttyACM0';
unix(command_check)
command_check = 'ls -l /dev/ttyACM1';
unix(command_check)

% Connect to spectroradiometer.
CMCheckInit(5); % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670


%% Measurement.
S=[380 2 201]; % Wavelength range (380-780 nm) / S=[380 2 201]
fw=MeasSpd(S,5,'all'); 
plot(SToWls(S),fw,'b-');


   % then measure 
    obj.radiometerObj.measure();

end