% %% 1) Spectrum Stability Test
% load('0721Spectrum(Stability).mat') % Load the data
% 
% % Wavelength linear interpolation (1/2/5 nm interval)
% t = [0:0.5:110.5]; % Time range in min
% w_2nm = [380:2:780]; % Given data
% w_1nm = [380:1:780]; % Interpolation range
% fw_1nm = interp1(w_2nm,fw,w_1nm);
% fw_5nm = fw_1nm(1:5:length(fw_1nm),:); % Data in 5 nm interval
% 
% load T_xyzJuddVos % Judd-Vos XYZ Color matching function
% T_XYZ = T_xyzJuddVos';
% XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation
% 
% time_min = 20; % Set the time range for comparison
% time_index = time_min*2;
% 
% figure(1); subplot(2,1,1); hold on;
% plot(w_2nm,fw(:,1:time_index),'k');
% plot(w_2nm,fw(:,time_index+1:222),'r');
% xlabel('Wavelength (nm)');
% ylabel('Spectral irradiance');
% % legend('');
% 
% figure(1); subplot(2,1,2); hold on;
% plot(t(1:time_index),XYZ(1:time_index,2),'k.-')
% plot(t(time_index+1:222),XYZ(time_index+1:222,2),'r.-')
% xlabel('Time (min)');
% ylabel('Luminance (cd/m2)');

%% 2) Spectrum Additivity Test
% This code compares the randomly generated spectra with the sum of single
% spectrum to see if its additivity works fine

% Initialize
clear all; close all; clc;

% Data preparation*********************************************************
date = 'Dataset8)';
add = '(WithinPrimary)_black10_Primary12_normal';

load(append(date,'Spectrum_single',add,'.mat')); % fw_single
load(append(date,'Spectrum_black',add,'.mat')); % fw_blk 
load(append(date,'Spectrum_white',add,'.mat')); % fw_white
load(append(date,'SpectrumInput',add,'.mat')) % SpectrumInput
load(append(date,'Spectrum_rand',add,'.mat')) % fw_rand

w = [380:2:780]; % Wavelength 380-780 nm with 2 nm interval

currentset = [252]; % Type the current level used
numcurrent = length(currentset); % num of measurements 
numsample = width(fw_rand); % Num of test samples

datatype = 1; % 1 = within primary / 2 = across primary
blkcorrect = 1; % 0 = raw spectra / 1 = black corrected

% Data should be typed till this part**************************************

% Single spectrum**********************************************************
% Black correction
if datatype == 1 ;
    
    % Black correction (WITHIN PRIMARY)
    fw_blks = ones(length(w),numcurrent); 
    fw_blks(:,:) = fw_blks.*fw_blk;

    if blkcorrect == 0  
            for i = 1:16;
            fw_singles{i} = fw_single(:,1+numcurrent*(i-1):numcurrent*(i))
            end
        elseif blkcorrect == 1
            for i = 1:16;
            fw_single_temp = fw_single(:,1+numcurrent*(i-1):numcurrent*(i))-fw_blks;
            fw_single_temp = max(fw_single_temp,0);
            fw_singles{i} = fw_single_temp;
            end
    end

elseif datatype == 2;
    
    % Black correction (ACROSS PRIMARY)
    fw_blks = ones(length(w),width(fw_single)); 
    fw_blks(:,:) = fw_blks.*fw_blk;

    if blkcorrect == 0  
       fw_single = fw_single;
    elseif blkcorrect ==1
       fw_single = fw_single - fw_blks;
       fw_single = max(fw_single,0); % Cut off the negative part
    end
    
    % Make cell with three separate components over primaries
    for p = 1:3
        fw_singles{p} = fw_single(:,1+16*(p-1):16*p);
    end     
end

% Calculate the sum of single spectra**************************************
if datatype == 1;
    
    % 1) Sum of single spectra (WITHIN PRIMARY)****************************
    % Find the step level number for the current (ex. current 252 = 23th out of 23)
    idx_SpectrumInput = zeros([size(SpectrumInput)]);
    current = [252];

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
    fw_blank = zeros(201,1);

    for s = 1:numsample % Num sample

        for ch=1:16 % Num channel

            fw_temp = cell2mat(fw_singles(ch)); 

            if idx_SpectrumInput(s,ch) == 0
               fw_rand_single(:,ch) = fw_blank;
            else
               fw_rand_single(:,ch) = fw_temp(:,idx_SpectrumInput(s,ch));        
            end

        end

        fw_rand_raw{:,s} = fw_rand_single;
        fw_sum = sum(fw_rand_single,2);
        fw_rand_sum(:,s) = fw_sum;
    end

elseif datatype == 2;
    
    % 2) Sum of single spectra (ACROSS PRIMARY)****************************
    % Find location of peak (LED channel 1-16)
    fw_R = cell2mat(fw_singles(1)); % Primary 1 / 201 x 16
    fw_G = cell2mat(fw_singles(2)); % Primary 2
    fw_B = cell2mat(fw_singles(3)); % Primary 3
    
    SpectrumInput_R = cell2mat(SpectrumInput(1)); % numsample x 16
    SpectrumInput_G = cell2mat(SpectrumInput(2));
    SpectrumInput_B = cell2mat(SpectrumInput(3));
    
    for i = 1:numsample
        idx_SpectrumInput_R(i,1) = find(SpectrumInput_R(i,:) == currentset); % 'numsample'x 1
        idx_SpectrumInput_G(i,1) = find(SpectrumInput_G(i,:) == currentset);
%         idx_SpectrumInput_B(i,1) = find(SpectrumInput_B(i,:) == currentset);
    end

    for s = 1:numsample % Num sample      

        fw_R_temp = zeros(201,16);
        fw_G_temp = zeros(201,16);
        fw_B_temp = zeros(201,16);
        
        fw_R_temp(:,idx_SpectrumInput_R(s)) = fw_R(:,idx_SpectrumInput_R(s));
        fw_G_temp(:,idx_SpectrumInput_G(s)) = fw_G(:,idx_SpectrumInput_G(s));
%         fw_B_temp(:,idx_SpectrumInput_B(s)) = fw_B(:,idx_SpectrumInput_B(s));

        fw_rand_R{s} = fw_R_temp;
        fw_rand_G{s} = fw_G_temp;
%         fw_rand_B{s} = fw_B_temp;Across
        
%         fw_sum = fw_R_temp + fw_G_temp + fw_B_temp; 
        fw_sum = fw_R_temp + fw_G_temp ; 
       
        fw_rand_raw{s} = fw_sum;
        fw_rand_sum(:,s) = sum(fw_sum,2);

    end
end

% Rand spectra********************************************************
% Black correction of the random mixed spectra
if blkcorrect == 0 % blackcorrection off
    fw_rand = fw_rand;
elseif blkcorrect == 1 % blackcorrection on
    fw_blks = ones(length(fw_rand),width(fw_rand));
    fw_blks = fw_blks.*fw_blk;
    fw_rand = fw_rand - fw_blks;
    fw_rand = max(fw_rand,0); % Cut off the negative parts
end

% Compare spectrum one-to-one**********************************************
numsample_part = 20;
start = 0; 

for i = 1:numsample_part
    figure(1); subplot(4,5,i); hold on;    
    plot(w,fw_rand(:,i+start),'k-')
    plot(w,fw_rand_sum(:,i+start),'r--')
    title(append('Test',num2str(i+start)));
    ylim([0 max(max([fw_rand fw_rand_sum]))]);
    % xlabel('Wavelength (nm)')
    % ylabel('Spectral irradiance')
end
legend('Measurement','Sum result','location','northeast')

%% Scale test 

figure(2); hold on;
plot(fw_rand,fw_rand_sum,'k.')
plot(fw_rand(:,2),fw_rand_sum(:,2),'r.','MarkerSize',8)
plot(fw_rand(:,3),fw_rand_sum(:,3),'r.','MarkerSize',8)
plot(fw_rand(:,6),fw_rand_sum(:,6),'r.','MarkerSize',8)
plot([0 max(max([fw_rand fw_rand_sum]))], [0 max(max([fw_rand fw_rand_sum]))],'k-')
xlabel('Radiance over wavelength (Single)')
ylabel('Radiance over wavelength (Sum)')


%% Compare one test color
testcolor = 2; % Test 2 3 6 (The others = 10)

fw_rand_raw_test = cell2mat(fw_rand_raw(testcolor));
fw_rand_R_test = cell2mat(fw_rand_R(testcolor));
fw_rand_G_test = cell2mat(fw_rand_G(testcolor));
% fw_rand_B_test = cell2mat(fw_rand_B(testcolor));
% testchannel = [idx_SpectrumInput_R(testcolor) idx_SpectrumInput_G(testcolor) idx_SpectrumInput_B(testcolor)];
testchannel = [idx_SpectrumInput_R(testcolor) idx_SpectrumInput_G(testcolor)];

fw_rand_test = fw_rand(:,testcolor);
fw_rand_sum_test = fw_rand_sum(:,testcolor);

ylimit = max(max([fw_rand_raw_test, fw_rand_test, fw_rand_sum_test]));

figure(10); subplot(4,1,1); hold on;
plot(w,fw_rand_raw_test);
plot(w,fw_rand_R_test,'r-');
plot(w,fw_rand_G_test,'g-');
% plot(w,fw_rand_B_test,'b-');
ylim([0 ylimit])
title('Raw spectra before summation')
legend(append('Channel ',num2str(testchannel)))

figure(10); subplot(4,1,2);
plot(w,fw_rand_sum_test,'r--');
ylim([0 ylimit])
title('Sum spectrum')

figure(10); subplot(4,1,3);
plot(w,fw_rand_test,'k-');
ylim([0 ylimit])
title('Measured spectrum')

figure(10); subplot(4,1,4); hold on;
plot(w,fw_rand_test,'k-');
plot(w,fw_rand_sum_test,'r--');
ylim([0 ylimit])
legend('Measure','Sum')
title('Sum vs. Measured spectrum')

%% XYZ Comparison
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range

fw_rand_1nm = interp1(w_2nm,fw_rand,w_1nm); % fw_rand
fw_rand_5nm = fw_rand_1nm(1:5:length(fw_rand_1nm),:); 

fw_rand_sum_1nm = interp1(w_2nm,fw_rand_sum,w_1nm); % fw_rand_sum
fw_rand_sum_5nm = fw_rand_sum_1nm(1:5:length(fw_rand_sum_1nm),:); 

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ_rand = 683*fw_rand_5nm'*T_XYZ; % XYZ calculation
XYZ_rand_sum = 683*fw_rand_sum_5nm'*T_XYZ;
xyY_rand = XYZToxyY(XYZ_rand');
xyY_rand_sum = XYZToxyY(XYZ_rand_sum');

% Colorgamut
colorgamut=XYZToxyY(T_XYZ');
colorgamut(:,82)=colorgamut(:,1);

% Spectrum
figure(4); subplot(2,2,1); hold on;
plot(w,fw_rand,'k-')
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
legend('Measurement (N=100)','location','northwest')
title('Spectrum (Measurement)')

figure(4); subplot(2,2,2); hold on;
plot(w,fw_rand_sum,'r-')
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
legend('Sum result (N=100)','location','northwest')
title('Spectrum (Sum result)')

% Luminance
numsampleset = [1:1:numsample];
figure(4);subplot(2,2,3); hold on;
plot(numsampleset,xyY_rand(3,:),'k.','markersize',7);
plot(numsampleset,xyY_rand_sum(3,:),'r+','markersize',7);
plot([numsampleset;numsampleset],[xyY_rand(3,:); xyY_rand_sum(3,:)],'k-'); 
xlabel('Test color')
ylabel('Luminance (cd/m2)')
title('Luminance')
legend('Measurement','Sum result','location','southeast');

% CIE (x,y) chromaticity
figure(4);subplot(2,2,4); hold on;
plot(xyY_rand(1,:),xyY_rand(2,:),'k.');
plot(xyY_rand_sum(1,:),xyY_rand_sum(2,:),'r+'); 
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
title('CIE (x,y) chromaticity')
legend('Measurement','Sum result');

%% Basic analysis on all channels
% Normalize the spectrum***************************************************
for ch=1:16
    fw_temp = cell2mat(fw_singles(ch));
    fw_norm(:,ch) = fw_temp./max(fw_temp);
end

% Test plot
% ch=1:16;
% figure(1); subplot(3,1,1); 
% plot(w,cell2mat(fw_singles(ch)))
% xlabel('Wavelength (nm)')
% ylabel('Spectral irradiance')
% title('Raw spectrum')
% 
% figure(1); subplot(3,1,2); 
% plot(w,fw_norm)
% xlabel('Wavelength (nm)')
% ylabel('Normalized Spectral irradiance')
% title('Normalized spectrum')
% 
% fw_temp = sum(fw_single,2)
% figure(1); subplot(3,1,3); hold on;
% plot(w,fw_white,'k-')
% plot(w,fw_temp,'r--')% Cut the negative parts off
% xlabel('Wavelength (nm)')
% ylabel('Spectral irradiance')
% title('White vs. Sum white')
% legend('white','Sum White')

% Peak spectrum************************************************************
for ch = 1:16; % CH = 1-16 (no data at 9)

    fw_temp = cell2mat(fw_singles(ch));
    % assignin('base','fw_temp',fw_single(:,1+repeat*(ch-1):repeat*(ch)))

    % Calculate the peak wavelength
    fw_peak = max(fw_temp);

    for i=1:numcurrent
    
        if sum(fw_temp(:,i))==0 
           continue % Pass the column with no measurement data 
        end
  
        idx_peak(i,ch) = find(fw_temp(:,i) == fw_peak(i));
        w_peak(i,ch) = 380+2*(idx_peak(i,ch)-1); % Starts from 380 nm with 2nm interval
    end
end 

% Find the half peak wavelength********************************************
for i=1:numcurrent

   for j=1:numcurrent
       if or(idx_peak(i) ==0, idx_peak(i) == 201) % Either no peak of noise
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

% XYZ calculation**********************************************************
w_2nm = [380:2:780]; % Given data
w_1nm = [380:1:780]; % Interpolation range
w_5nm = [380:5:780];
fw_1nm = interp1(w_2nm,fw_temp,w_1nm);
fw_1nm = fw_1nm';
fw_5nm = fw_1nm(1:5:length(w_1nm),:); % Data in 5 nm interval

load T_xyzJuddVos % Judd-Vos XYZ Color matching function
T_XYZ = T_xyzJuddVos';
XYZ = 683*fw_5nm'*T_XYZ; % XYZ calculation
xyY = XYZToxyY(XYZ');
colorgamut=XYZToxyY(T_XYZ');X = [1 2 3 4; 5 6 7 8; 9 10 11 12];
colorgamut(:,82)=colorgamut(:,1);

figure(2); subplot(2,2,1); hold on;
plot(w,fw_temp)
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
title('Spectrum')

% Add legend with increasing number
for i=1:length(currentset)
  leg{i}=sprintf('current %d',i);
end
% legend(leg,'fontsize',6);

% CIE (x,y) chromaticity
figure(2);subplot(2,2,2); hold on;
plot(xyY(1,:),xyY(2,:),'r.'); % Measurement point
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
title('CIE (x,y) chromaticity')
legend('Test');

% Wavelength
figure(2);subplot(2,2,3); hold on;
plot(currentset,w_peak,'g.');
xlabel('Current')
ylabel('Wavelength (nm)')
xlim([0 252]);
ylim([380 780]);
title('Peak Wavelength')

% Bandwidth
figure(2);subplot(2,2,4); hold on;X = [1 2 3 4; 5 6 7 8; 9 10 11 12];
plot(currentset,w_bandwidth,'b.');
xlabel('Current')
ylabel('Bandwidth (nm)')
xlim([0 252]);
ylim([0 max(w_bandwidth)+5])
title('Bandwidth')

%% Error over the SpectrumInput (For Within Primary comparison only for this version)
% fw_error = abs(fw_rand_sum - fw_rand);
% error = sum(fw_error);
% 
% findvalue = 0;
% for i=1:numsample
%     numfindvalue = sum(SpectrumInput(i,:)==findvalue);
%     numfindvalueset_0(i) = numfindvalue;
% end
% 
% findvalue = 252;
% for i=1:numsample
%     numfindvalue = sum(SpectrumInput(i,:)==findvalue);
%     numfindvalueset_252(i) = numfindvalue;
% end
% 
% sumSpectrumInput = sum(SpectrumInput,2);
% 
% figure(7); subplot(2,2,1); hold on; 
% plot(numfindvalueset_0,error,'b.')
% xlabel('Number of LED channel')
% ylabel('Error')
% legend('current=0')
% title('Error vs. Num channel-1')
% 
% figure(7); subplot(2,2,2); hold on; 
% plot(numfindvalueset_252,error,'b.')
% xlabel('Number of LED channel')
% ylabel('Error')
% legend('current=252')
% title('Error vs. Num channel-2')
% 
% figure(7); subplot(2,2,3); hold on;
% plot(sumSpectrumInput,error,'r.')
% xlabel('Sum channel currents')
% ylabel('Error')
% legend('(N=100)')
% title('Error vs. Sum channel currents')