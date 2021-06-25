%% Camera MTF measurement 
% Written by Semin Oh
% As of June 11th

% Initialize
clear; close all; clc;

% DMD rectangle size (pre-calculated)
boundingbox_half_right = [1564.5,953.5,195,219]; % top left (x,y) coordinates of the half right side of the test pattern  [x, y, horizontal, vertical]
boudingbox = [1564.5-195,953.5,195*2,219]; % DMD size [x, y, width, height]

% Load image and crop the target image 
image = imread('blkandwhite.tiff','tiff');
[Ypixel Xpixel] = size(image);

imagecrop = image(boudingbox(2):boudingbox(2)+boudingbox(4),boudingbox(1):boudingbox(1)+boudingbox(3));
[Ypixel_crop Xpixel_crop] = size(imagecrop);

% Resize the cropped image as desired (clear target)
imagecrop = imagecrop(0.1*Ypixel_crop:0.9*Ypixel_crop,0.1*Xpixel_crop:0.9*Xpixel_crop);
[Ypixel_crop Xpixel_crop] = size(imagecrop);

% figure(1); imshow(image); title('Original Image');
figure(2); imshow(imagecrop); title('Cropped Image');


% Draw 25% / 50% / 75% position of the cropped image
imagecrop_25 = imagecrop(round(0.25*Ypixel_crop),:);
imagecrop_50 = imagecrop(round(0.50*Ypixel_crop),:);
imagecrop_75 = imagecrop(round(0.75*Ypixel_crop),:);
imagecrop_avg = mean([imagecrop_25;imagecrop_50;imagecrop_75]);
% imagecrop_std = std([imagecrop_25;imagecrop_50;imagecrop_75]);

white = max(imagecrop_avg);
black = min(imagecrop_avg);
contrast = (white-black)/(white+black)

figure(3); hold on;
plot(1:Xpixel_crop,imagecrop_25(:),'k-');
plot(1:Xpixel_crop,imagecrop_50(:),'g-');
plot(1:Xpixel_crop,imagecrop_75(:),'b-');
plot(1:Xpixel_crop,imagecrop_avg(:),'r--','LineWidth',2);

% error-bar for the average grap
% errorindex = 20;
% x = 1:errorindex:Xpixel_crop;
% y = imagecrop_all(1:errorindex:Xpixel_crop);
% err = imagecrop_std(1:errorindex:Xpixel_crop);
% errorbar(x,y,err,'horizontal','o')
title('Line spread function');
xlabel('Pixels (Horizontal)');
ylabel('Intensity');
legend('25%','50%','75%','Avg');

txt1 = num2str(contrast);
text(350,50,txt1,'FontSize',14,'Color','r','HorizontalAlignment','right')


%% 1) MTF measure
% Calculate MTF (modulus of OTF)
FFT_25 = abs(fftshift(fft(imagecrop_25)));
FFT_50 = abs(fftshift(fft(imagecrop_50)));
FFT_75 = abs(fftshift(fft(imagecrop_75)));
FFT_avg = abs(fftshift(fft(imagecrop_avg)));

% Normalize to max value as 1
FFT_25 = FFT_25/max(FFT_25);
FFT_50 = FFT_50/max(FFT_50);
FFT_75 = FFT_75/max(FFT_75);
FFT_avg = FFT_avg/max(FFT_avg);

% Set x-axis (frequecy) *needs to be updated
Size = length(FFT_25);
spacing = 20;% spacing between data points 
fsx = abs(1/spacing); % turn into sampling frequency
a = linspace(-Size/2,Size/2,Size); % form scale for conversion based on frequency bins
conversionx = fsx./Size; % conversion factor for frequency bin units to frequency (unit^-1)
Psi = a.*conversionx; % frequency (unit^-1)

% Plot
figure(4); hold on;
plot(Psi,FFT_25(:),'k-');
plot(Psi,FFT_50(:),'g-');
plot(Psi,FFT_75(:),'b-');
plot(Psi,FFT_avg(:),'r--');
xlabel('Frequency distribution');
ylabel('');
legend('25%','50%','75%','Avg');

sumFFT_25 = sum(FFT_25);
sumFFT_50 = sum(FFT_50);
sumFFT_75 = sum(FFT_75);
sumFFT_avg = sum(FFT_avg)

txt2 = num2str(sumFFT_avg);
text(0.02,0.5,txt2,'FontSize',14,'Color','r','HorizontalAlignment','right')

%% correct sampling frequency for conversion on frequency bins to frequency
% figure(5); hold;
% Size = length(MTF_25);
% spacing = 100;% spacing between data points 
% fsx = abs(1/spacing); % turn into sampling frequency
% a = linspace(-Size/2,Size/2,Size); % form scale for conversion based on frequency bins
% conversionx = fsx./Size; % conversion factor for frequency bin units to frequency (unit^-1)
% Psi = a.*conversionx; % frequency (unit^-1)
% plot(Psi,MTF_25)
% xlim([0,Psi(end)])
% xlabel('cycles/unit')
% ylabel('MTF')
