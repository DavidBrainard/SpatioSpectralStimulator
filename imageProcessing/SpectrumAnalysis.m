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

%% 2) Spectrum Additivity
%% Load Single spectrum
clear; close all; clc;

load('0722Spectrum_single.mat'); % fw_single
blacklevel = 10;
load(append('black',num2str(blacklevel),'(ch6-15)','.mat')) % fw 
fw_blk = fw;

w = [380:2:780];
current = [30:10:240,252]; % measured currents level 30-240 with 10 interval, 252 (N=23)
numcurrent = 23; % num of measurements 

% Save all single channels in a cell (Ch.1-16 / VPixx Ch.0-15)
fw_blks_rand = ones(length(w),numcurrent); 
fw_blks_rand(:,:) = fw_blks_rand.*fw_blk;

blkcorrect = 1 % 0 = raw spectra / 1 = black corrected

if blkcorrect == 0  
        for i = 1:16;
        fw_singles{i} = fw_single(:,1+numcurrent*(i-1):numcurrent*(i))
        end
    elseif blkcorrect == 1
        for i = 1:16;
        fw_singles{i} = fw_single(:,1+numcurrent*(i-1):numcurrent*(i))-fw_blks_rand;
        end
end

% Cut the negative parts off
for i=1:16
    temp = cell2mat(fw_singles(i)); % 'temp' is 201x23 double mat
    for j=1:numcurrent
        for k=1:201
            if temp(k,j)<0
                temp(k,j)=0;
            end
        end
    end
    fw_singles{i} = temp
end

% Test plot
ch=16;
plot(w,cell2mat(fw_singles(ch)))
%% Peak and half-peak calculations (by channel)

ch=10% 1-16 (no data at 9)
%CH 1-5 are missing / CH14-16 = 30 / CH 12-13 = 40 / CH11 = 60 / CH10 = 40 / CH8 = 70 / CH7 =120 / CH6 = 40

fw_temp = cell2mat(fw_singles(ch));
% assignin('base','fw_temp',fw_single(:,1+repeat*(ch-1):repeat*(ch)))

% Calculate the peak wavelength
fw_peak = max(fw_temp);

for i=1:numcurrent
    
    if sum(fw_temp(:,i))==0 
       continue % Pass the no measurement data column
    end
  
    idx_peak(i) = find(fw_temp(:,i) == fw_peak(i));
    w_peak(i) = 380+2*(idx_peak(i)-1); % Starts from 380 nm with 2nm interval
    
end

% Find the half peak wavelength
fw_norm = fw_temp./max(fw_temp)

for i=1:numcurrent

   for j=1:numcurrent
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
fw_1nm = interp1(w_2nm,fw_temp,w_1nm);
fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation
xyY = XYZToxyY(XYZ');
colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

figure(1); subplot(2,2,1); hold on;
plot(w,fw_temp)
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

%% Load rand spectra 
load('0723SpectrumInput.mat') % SpectrumInput
load('0723Spectrum_rand.mat') % fw_rand

w = [380:2:780];
numsample = width(fw_rand);

% Black correction of the random mixed spectra
blkcorrect = 1; % 0=off / 1=on

if blkcorrect == 0 % blackcorrection off
    fw_rand = fw_rand;
elseif blkcorrect == 1 % blackcorrection on
    fw_blks_rand = ones(length(fw_rand),width(fw_rand));
    fw_blks_rand = fw_blks_rand.*fw_blk;
    fw_rand = fw_rand - fw_blks_rand;
end

% Cut the negative parts off
for j=1:numsample
   for k=1:201
      if fw_rand(k,j)<0
         fw_rand(k,j)=0;
      end
   end
end
    
%% Sum of single spectra
% Random spectra was made using the current levels among four off/150/200/252
% Current level = 30-240 in 10 interval, and 252 as maximum

% Find the step level number for the current (ex. current 252 = 23th out of 23)
idx_SpectrumInput = zeros([size(SpectrumInput)]);

for j=1:numsample % Num samples
    for i = 1:16 % Num channels
        idx_current = find(current==SpectrumInput(j,i));
        
            % Ignore the empty cells (replace it as 0)
            if isempty(idx_current) == 1;
                idx_current = 0;
            else 
                idx_current = idx_current;
            end
            
        idx_SpectrumInput(j,i) = idx_current;
    end
end

% Load current data from each single peak
fw_rand_single = zeros(201,16);

for j=1:numsample % Num sample
    for i=1:16 % Num channel
       
         if idx_SpectrumInput(j,i) == 0; 
           continue % Pass the remaining command
         end
                                
        fw_temp = cell2mat(fw_singles(i)); % 201x23 mat
        fw_rand_single(:,i) = fw_temp(:,idx_SpectrumInput(j,i));        
     end
    
    fw_sum = sum(fw_rand_single,2);
    fw_rand_sum(:,j) = fw_sum;
end

%% Plot the graph (current = off/150/200/252)
figure(1); subplot(2,1,1); hold on;

ch=8;

plot(SpectrumInput(:,ch),'r-')
plot(1:16,SpectrumInput,'+');
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

%% Comparison
figure(2); subplot(2,1,1); hold on;
plot(w,fw_rand,'k-')
plot(w,fw_rand_sum,'r--')
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')

sample=6;
figure(2); subplot(2,1,2); hold on;
plot(w,fw_rand(:,sample),'k-')
plot(w,fw_rand_sum(:,sample),'r--')
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
