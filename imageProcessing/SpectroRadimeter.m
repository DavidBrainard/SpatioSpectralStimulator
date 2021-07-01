%% VPixx - setting the ground state of the projector

% Initialize
clear all; close all; clc;
tbUse('BrainardLabBase')

% Run terminal command 'vputil' (case sensitive) / both unix and system
% work
command = 'vputil';
unix(command)

devsel ppx % Select the PROPixx with the command (Just type this on the command)

% VPixx Primary settings (type each on the command screen)
rw 0x1c8 0x1 % (Primary 0 always ON) > Screen goes red 
rw 0x1c8 0x2 % (Primary 1 always ON) > Screen goes green
rw 0x1c8 0x4 % (Primary 2 always ON) > Screen goes blue
rw 0x1c8 0x7 % (All primaries on) > Screen goes brighter with no chromaticity changes (at least visually)
rw 0x1c8 0x0 % (Set to default) > Back to the first screen

%% Measurement using PR-670 Spectroradiometer
% In terminal, command to authorize to read and write 'sudo chomd a+rw /dev/ttyACM0'
% Also command to check the status 'ls -l /dev/ttyACM0'

CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670
S=[380 2 201];
fw=MeasSpd(S,5,'all'); % Measurement (it takes a while)
plot(SToWls(S),fw);
xlabel('Wavelength(nm)')
ylabel('Spectral irradiance')

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
load T_cones_sp % Smith-Pokorny Cone spectral sensitivity function
T_Cones = T_cones_sp';

XYZ = fw'*T_XYZ; % XYZ calculation
Cones = fw'*T_Cones; % Cone signals calculation
xyY = XYZToxyY(XYZ); % xyY
XYZw = [100 100 100]'; % Lab white
Lab = XYZToLab(XYZ,XYZw); % Lab calculation

