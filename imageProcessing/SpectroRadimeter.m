%% VPixx - setting the ground state of the projector

% Initialize
clear all; close all; clc;
tbUse('BrainardLabBase')

%% Measurement using PR-670 Spectroradiometer
% In terminal, command to authorize to read and write 'sudo chmod a+rw /dev/ttyACM0 (the last character is number zero)'
% Also command to check the status 'ls -l /dev/ttyACM0'

%% Initialize (authorization)
% Mostly PR670 connects to 'ttyACM0', but sometimes to 'ttyACM1', so just
% do both to make it sure to be connected

command_PR670 = 'sudo chmod a+rw /dev/ttyACM0';
unix(command_PR670)

command_PR670 = 'sudo chmod a+rw /dev/ttyACM1';
unix(command_PR670)

command_check = 'ls -l /dev/ttyACM0';
unix(command_check)

%% Measurement 
CMCheckInit(5) % CMCheckInit([meterType], [PortString]) / 5 is allocated to PR670

S=[380 2 201]; % Wavelength range (380-780 nm) / S=[380 2 201]
fw=MeasSpd(S,5,'all'); % Measurement (it takes a while)
plot(SToWls(S),fw,'b-');

num=subColor;
filename = append('sub',num2str(num),'.mat'); 
% filename = append('blk','.mat'); 
save(filename,'fw')

%%
load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw'*T_XYZ; % XYZ calculation
xyY = XYZToxyY(XYZ'); % xyY
colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

XYZw = [100 100 100]'; % Lab white
Lab = XYZToLab(XYZ',XYZw); % Lab calculation

% Power spectrum
figure(1); subplot(2,2,1); hold on;
plot(SToWls(S),fw,'b-');
xlabel('Wavelength(nm)')
ylabel('Spectral irradiance')
xlim([380 780]);

% CIE (x,y) chromaticity
figure(1);subplot(2,2,2); hold on;
plot(xyY(1,:),xyY(2,:),'r.'); % Measurement point
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
legend('Test');

% CIELAB 
figure(1);subplot(2,2,3); hold on;
plot(Lab(1,:),Lab(2,:),'r.'); % Measurement point
xlabel('CIELAB a*')
ylabel('CIELAB b*')
xlim([-100 100]);
ylim([-100 100]);
legend('Test');


%% Cone signals calculation
% load T_cones_sp % Smith-Pokorny Cone spectral sensitivity function
% T_Cones = T_cones_sp';
% Cones = fw'.*T_Cones; % Cone signals calculation

