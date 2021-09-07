%% Load variables with the same name (variable name = 'fw')

% Initialize
clear all; close all; clc;

% Set wavelength range 380-780 nm with 5 nm interval
w = [380:5:780];

% Load white and black variables
white = load(append('white','.mat'));
black = load(append('black','.mat'));
white_no = load(append('white','(no)','.mat'));
black_no = load(append('black','(no)','.mat'));

% Read the spectrum
spd_w = white.white;
spd_blk = black.black;
spd_w_no = white_no.fw;
spd_blk_no = black_no.fw;

% Normalize its max value to 1
spd_w = spd_w./max(spd_w);
spd_blk = spd_blk./max(spd_blk);
spd_w_no = spd_w_no./max(spd_w_no);
spd_blk_no = spd_blk_no./max(spd_blk_no);

% LED Spectrums with the pinhole (channel 3 and 8 are missing / cannot
% measure)
for i=1:16
    subcolor(i) = load(append('subcolor_',num2str(i-1),'.mat'));
    subcolors(:,i) = subcolor(i).fw;
    subcolors_norm(:,i) = subcolors(:,i)./max(subcolors(:,i));
end

% LED Spectrums without the pinhole
for i=1:16
    subcolor(i) = load(append('subcolor_',num2str(i-1),'(no)','.mat'));
    subcolors_no(:,i) = subcolor(i).fw;
    subcolors_no_norm(:,i) = subcolors_no(:,i)./max(subcolors_no(:,i));
end

% Plot the spectrum
figure(1); hold on;subplot(3,1,1);
plot(w,subcolors_norm);
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
ylim([0 1]);
title('LED 16-channels with Pinhole (5.8mm)');

figure(1); hold on;subplot(3,1,2);
plot(w,subcolors_no_norm);
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
ylim([0 1]);
title('LED 16-channels without Pinhole');

figure(1); subplot(3,1,3); hold on; 
plot(w,spd_w,'k-'); 
plot(w,spd_w_no,'g--');
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
ylim([0 1]);
title('White (all 16-channels on)');
legend('Pinhole','No pinhole');

%% XYZ calculation 
load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*subcolors'*T_XYZ; % XYZ calculation
XYZ_no = 683*subcolors_no'*T_XYZ;
xyY = XYZToxyY(XYZ'); % xyY
xyY_no = XYZToxyY(XYZ_no');

colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

% CIELAB calculation
Y = max(xyY(3,:));
XYZw = [Y Y Y]'; % Lab white (Arbitrary set to equi-energy)
Y_no = max(xyY_no(3,:));
XYZw_no = [Y_no Y_no Y_no]';
Lab = XYZToLab(XYZ',XYZw); % Lab calculation
Lab_no = XYZToLab(XYZ_no',XYZw_no);

% Power spectrum
figure(2); subplot(2,2,1:2); hold on;
plot(w,subcolors_norm,'k-');
plot(w,subcolors_no_norm,'g--');
xlabel('Wavelength(nm)')
ylabel('Spectral irradiance')
xlim([380 780]);
legend('Pinhole','No pinhole');
title('Spectrum');

% CIE (x,y) chromaticity
figure(2);subplot(2,2,3); hold on;
plot(xyY(1,:),xyY(2,:),'k.'); % Measurement point
plot(xyY_no(1,:),xyY_no(2,:),'g.'); % Measurement point
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
legend('Pinhole','No pinhole');
title('CIE (x, y) chromaticity');

% CIELAB
figure(2);subplot(2,2,4); hold on; grid on;
plot(Lab(2,:),Lab(3,:),'k.');
plot(Lab_no(2,:),Lab_no(3,:),'g.');
xlabel('CIELAB a*')
ylabel('CIELAB b*')
xlim([-200 200]);
ylim([-200 200]);
legend('Pinhole','No pinhole','location','southeast');
title('CIELAB (a*, b*) plane');