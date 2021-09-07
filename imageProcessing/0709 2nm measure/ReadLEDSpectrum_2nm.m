%% Load variables with the same name (variable name = 'fw')
% This file is for reading the LED spectrum files in 2nm interval
% as of 0709 measurements

% Initialize
clear all; close all; clc;

% Set wavelength range 380-780 nm with 5 nm interval
w = [380:2:780];

% Load white and black variables
white = load(append('white','.mat'));
black = load(append('black','.mat'));

% Read the spectrum
spd_w = white.fw;
spd_blk = black.fw;

% Normalize its max value to 1
spd_w_norm = spd_w./max(spd_w);
spd_blk_norm = spd_blk./max(spd_blk);

% LED Spectrums with the pinhole (channel 3 and 8 are missing / cannot
% measure)
for i=1:16
    subcolor(i) = load(append('sub',num2str(i-1),'.mat'));
    spd_raw(:,i) = subcolor(i).fw;
    spd_norm(:,i) = spd_raw(:,i)./max(spd_raw(:,i));
end

% Erase the nineth channel which doesn't contain any data
peak_spd = max(spd_raw);
peak_spd(:,9) = [];
spd_raw(:,9) = []; 
spd_norm(:,9) = [];

spd_rgbw = sum(spd_raw,2); % White - Sum of all 15 channels

% Calculate the peak wavelength
for i=1:15
    idx_peak(i) = find(spd_raw(:,i) == peak_spd(i));
    w_peak(i) = 380+2*(idx_peak(i)-1); % Starts from 380 nm with 2nm interval
end

% Find the half peak wavelength
% [val,idx] = min(abs(spd_norm-0.5));
% w_halfpeak = 380+2*(idx-1);
% w_diff_halfpeak = abs(w_peak - w_halfpeak);
% 
% w_halfpeak_L = w_peak - w_diff_halfpeak;
% w_halfpeak_R = w_peak + w_diff_halfpeak;

% Find the ground touching wavelength (20% peak)
% [val2,idx2] = min(abs(spd_norm-0.2));
% w_ground = 380+2*(idx2-1);
% w_diff_ground = abs(w_peak - w_ground);
% 
% w_ground_L = w_peak - w_diff_ground;
% w_ground_R = w_peak + w_diff_ground;

% LED single spectrum power meter measure (as of 0715)
peak_power_nW = [32.775 45.885 85.250 75.492 173.250 480.753 320.259 525.146 596.578 279.033 678.067 767.097 1151 1368 912.877]; % Unit: nano watt
peak_power_mW = peak_power_nW.*0.001; % Unit: micro watt

peak_power_mW_norm = (peak_power_mW./max(peak_power_mW)).*max(peak_spd);

% Create the LED spectrum based on the power measure (Newport 1936-C)
spd_power = zeros(201,15);
for i=1:15
    spd_power(idx_peak(i),i) = peak_power_mW(i);
end

spd_power_norm = (spd_power./max(max(spd_power))).*max(max(spd_raw));

% Normalized LED spectrum
figure(1); hold on;subplot(3,1,1);
plot(w,spd_norm);
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
ylim([0 1]);
title('LED 16-channels with Pinhole (5.8mm)');

% Raw spectrum
figure(1); subplot(3,1,2); hold on; 
plot(w,spd_raw);
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
title('LED 16-channels with Pinhole (5.8mm)');

% White spectrum
figure(1); subplot(3,1,3); hold on; 
plot(w,spd_w_norm,'k-'); 
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
xlim([380 780]);
ylim([0 1]);
title('White (all 16-channels on)');
legend('white');

% Calculation of dependency error
error_dependency = 100*(((peak_spd)-(peak_power_mW_norm))./((peak_spd+peak_power_mW_norm)./2));

%% LED Spectrum devince dependency (PR670 vs. Powermeter 1936-C)

% Overview
figure(2); hold on;
plot(w,spd_raw,'k-');
plot(w,spd_power_norm,'r-');
plot(w_peak,peak_power_mW_norm,'ro','MarkerSize',4,'MarkerFaceColor','k');
xlim([380 780]);
ylim([0 max(max(spd_raw))]);
xlabel('Wavelength (nm)','fontsize',14)
ylabel('Spectral irradiacne','fontsize',14);
title('LED Spectrum Measurement device dependency');
legend('PR670','','','','','','','','','','','','','','','','1936-C','location','northwest','fontsize',13);

% View by each channel
figure(3); hold on;
for i=1:15
    subplot(3,5,i); hold on;
    plot(w,spd_raw(:,i),'k-');
    plot(w,spd_power_norm(:,i),'r-');
    xlim([380 780]);
    ylim([0 max(max(spd_raw))]);
    title(append('Ch',num2str(i),' (',num2str(w_peak(i)),' nm)'));
end
xlabel('Wavelength (nm)');
ylabel('Spectral irradiacne');
legend('PR670','1936-C','location','northwest');


%% LED Spectrum additivity (white vs. sum white)

figure(4); hold on; 
plot(w,spd_w,'k-','linewidth',1);
plot(w,spd_rgbw,'g--','linewidth',1);
xlabel('Wavelength (nm)','fontsize',14);
ylabel('Spectral irradiacne (raw)','fontsize',14);
xlim([380 780]);
title('LED Spectrum additivity');
legend('white','Sum of all 15 channels','location','northwest','fontsize',13);

error_additivity = 100*(abs(spd_w - spd_rgbw)./spd_w)

% figure(5);
% plot(w,error_additivity);

%% XYZ calculation 
% load T_xyzJuddVos % Judd-Vos XYZ Color matching function
% T_XYZ = T_xyzJuddVos';
% XYZ = 683*subcolors'*T_XYZ; % XYZ calculation
% XYZ_no = 683*subcolors_no'*T_XYZ;
% xyY = XYZToxyY(XYZ'); % xyY
% xyY_no = XYZToxyY(XYZ_no');
% 
% colorgamut=XYZToxyY(T_XYZ');
% colorgamut(:,82)=colorgamut(:,1);
% 
% % CIELAB calculation
% Y = max(xyY(3,:));
% XYZw = [Y Y Y]'; % Lab white (Arbitrary set to equi-energy)
% Y_no = max(xyY_no(3,:));
% XYZw_no = [Y_no Y_no Y_no]';
% Lab = XYZToLab(XYZ',XYZw); % Lab calculation
% Lab_no = XYZToLab(XYZ_no',XYZw_no);
% 
% % Power spectrum
% figure(2); subplot(2,2,1:2); hold on;
% plot(w,subcolors_norm,'k-');
% plot(w,subcolors_no_norm,'g--');
% xlabel('Wavelength(nm)')
% ylabel('Spectral irradiance')
% xlim([380 780]);
% legend('Pinhole','No pinhole');
% title('Spectrum');
% 
% % CIE (x,y) chromaticity
% figure(2);subplot(2,2,3); hold on;
% plot(xyY(1,:),xyY(2,:),'k.'); % Measurement point
% plot(xyY_no(1,:),xyY_no(2,:),'g.'); % Measurement point
% plot(colorgamut(1,:),colorgamut(2,:),'k-');
% xlabel('CIE x')
% ylabel('CIE y')
% xlim([0 1]);
% ylim([0 1]);
% legend('Pinhole','No pinhole');
% title('CIE (x, y) chromaticity');
% 
% % CIELAB
% figure(2);subplot(2,2,4); hold on; grid on;
% plot(Lab(2,:),Lab(3,:),'k.');
% plot(Lab_no(2,:),Lab_no(3,:),'g.');
% xlabel('CIELAB a*')
% ylabel('CIELAB b*')
% xlim([-200 200]);
% ylim([-200 200]);
% legend('Pinhole','No pinhole','location','southeast');
% title('CIELAB (a*, b*) plane');