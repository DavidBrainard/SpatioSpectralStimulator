function [] = SetProjectorSubprimarySettings(testSubprimarySettings,options)
% Set the projector subprimary settings for SACC project.
%
% Syntax: [] = SetProjectorSubprimarySettings(testSubprimarySettings)
%
% Description:
%    This sets up the projector subprimary settings for the project. This
%    should be used wherever the projector displays an image including
%    calibration, displaying gabor patch, etc.
%
% Inputs:
%    testSubprimarySettings -      Projector input settings that reproduce
%                                 the desired contrast.
%
% Optional key/value pairs:
%    'projectorMode' -            Boolean (default true). Set the projector
%                                 pulse mode either 'Normal' (true) or
%                                 'Steady-on' (false).
%    'measurementOption' -        Boolean (default true). Set if you want
%                                 to proceed the measurements. If you set
%                                 'true', measurement will be included, and
%                                 'false' will skip the measurement. This
%                                 will be useful if you run the code
%                                 outside the lab and debugging the code.
%    'verbose' -                  Boolean. Default true.  Controls plotting
%                                 and printout.
%
% History:
%    10/28/21  smo                Started on it

%% Set parameters.
arguments
    testSubprimarySettings
    options.projectorMode (1,1) = true
    options.measurementOption (1,1) = true
    options.verbose (1,1) = true
end

%% Connect the projector and set the projector mode.
% Add VPixx toolbox 'Datapixx' to the path.
toolboxDirectory = '/home/colorlab/Documents/MATLAB/toolboxes/VPixx'; % This is where Linux box store the file.
addpath(genpath(toolboxDirectory));

% Connect to the Vpixx projector.
isReady = Datapixx('open');
isReady = Datapixx('IsReady');

% Set the projector mode.
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