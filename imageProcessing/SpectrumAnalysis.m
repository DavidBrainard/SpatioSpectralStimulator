%% 1) Stability over time 
% Save the XYZ values
% Wavelength linear interpolation (1/2/5 nm interval)
load('0721Spectrum(Stability).mat') % Load the data

t = [0:0.5:110.5]; % Time range in min
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range
fw_1nm = interp1(w_2nm,fw,w_1nm);
fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation

time_min = 20; % Set the time range for comparison
time_index = time_min*2;

figure(1); subplot(2,1,1); hold on;
plot(w_2nm,fw(:,1:time_index),'k');
plot(w_2nm,fw(:,time_index+1:222),'r');
xlabel('Wavelength (nm)');
ylabel('Spectral irradiance');
% legend('');

figure(1); subplot(2,1,2); hold on;
plot(t(1:time_index),XYZ(1:time_index,2),'k.-')
plot(t(time_index+1:222),XYZ(time_index+1:222,2),'r.-')
xlabel('Time (min)');
ylabel('Luminance (cd/m2)');

%% 2) Single spectrum analysis
clear; close all; clc;

load('0722Spectrum_single.mat');
w = [380:2:780];
current = [30:10:240,252];

% plot(w,fw_single)

ch=11% 1-16 (no data at 9)
%CH 1-5 are missing / CH14-16 = 30 / CH 12-13 = 40 / CH11 = 60 / CH10 = 40 / CH8 = 70 / CH7 =120 / CH6 = 40

assignin('base','fw_onechannel',fw_single(:,1+23*(ch-1):23*(ch)))
fw_onechannel;

% Calculate the peak wavelength
fw_peak = max(fw_onechannel);

for i=1:23
    
  for j=1:23
    if sum(fw_onechannel(:,i))==0 % Pass the no measurement data column
        i=i+1;
    end
  end
    
    idx_peak(i) = find(fw_onechannel(:,i) == fw_peak(i));
    w_peak(i) = 380+2*(idx_peak(i)-1); % Starts from 380 nm with 2nm interval
    
end

% Find the half peak wavelength
fw_norm = fw_onechannel./max(fw_onechannel)

for i=1:23

   for j=1:23
    if idx_peak(i) ==0
        i = i+1;
    end
   end
    
[val_L,idx_L] = min(abs(fw_norm(1:idx_peak(i),:)-0.5));
[val_R,idx_R] = min(abs(fw_norm(idx_peak(i):201,:)-0.5));

idx_R_true(i) = idx_peak(i) + idx_R(i);

end

w_halfpeak_L = 380+2*(idx_L-1);
w_halfpeak_R = 380+2*(idx_R_true-1);

w_bandwidth = w_halfpeak_R - w_halfpeak_L;

% XYZ calculation
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range
fw_1nm = interp1(w_2nm,fw_onechannel,w_1nm);
fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation
xyY = XYZToxyY(XYZ');
colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

figure(1); subplot(2,2,1); hold on;
plot(w,fw_onechannel)
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
title('Spectrum')

s = [30:10:240,252];
for i=1:length(s)
  leg{i}=sprintf('current %d',i);
end
% legend(leg,'fontsize',6);

% CIE (x,y) chromaticity
figure(1);subplot(2,2,2); hold on;
plot(xyY(1,:),xyY(2,:),'r.'); % Measurement point
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
title('CIE (x,y) chromaticity')
legend('Test');

% Wavelength
figure(1);subplot(2,2,3); hold on;
plot(current,w_peak,'g.');
xlabel('Current')
ylabel('Wavelength (nm)')
xlim([0 252]);
ylim([380 780]);
title('Peak Wavelength')

% Bandwidth
figure(1);subplot(2,2,4); hold on;
plot(current,w_bandwidth,'b.');
xlabel('Current')
ylabel('Bandwidth (nm)')
xlim([0 252]);
ylim([0 max(w_bandwidth)+5])
title('Bandwidth')

%% 3) Additivity test 
load('0723SpectrumInput.mat') % SpectrumInput
load('0723Spectrum_rand.mat') % fw_rand

w = [380:2:780];

% Spectra input (current = off/150/252)
figure(1); subplot(2,1,1); hold on;
SpectrumInput(:,9)=[];

ch=8;

plot(SpectrumInput(:,ch),'r-')
plot(1:15,SpectrumInput,'+');
xlabel('LED Channel')
ylabel('LED Current')
xlim([6 15])
legend('Test')

% Random spectra (N=30)
figure(1); subplot(2,1,2); hold on;
plot(w,fw_rand,'k-')
plot(w,fw_rand(:,ch),'r','LineWidth',2)
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
