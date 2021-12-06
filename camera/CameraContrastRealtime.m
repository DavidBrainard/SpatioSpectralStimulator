%% 1) Camera preview
clear; close all; clc;

% Load the camera info
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview
windowsize = 0.25; % The screen ratio to original. 0.25 means 25% raw image (type 0 to 1)
vidRes_resize = vidRes.*windowsize;

hFig = figure('Units', 'pixels', 'Position', [100 100 vidRes_resize(1) vidRes_resize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidRes_resize(1) vidRes_resize(2)]);
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

% Info about the marker position
imWidth = vidRes(1);
imHeight = vidRes(2);
numBands = vid.NumberOfBands;
markerindex = 0.5; % Set the centered point

% Camera preview with the centered point marked
figure(1);
preview(vid, hImage)
hLine = line(hAxes, round([markerindex*imWidth, markerindex*imWidth]),round([markerindex*imHeight, markerindex*imHeight]),'Marker','+','MarkerSize',30,'color','r','LineWidth',1,'LineStyle','none');
a = 0.45*imWidth; % rect starting x-coordinate
b = 0.46*imHeight; % rect starting y-coordinate / It is coefficient of 'a' variable in the following script
c = 309; % rect width
d = 166; % rect height
rectangle('Position',[a,b,c,d],'Curvature',[0,0],'LineWidth',1,'LineStyle','--','edgecolor','y')
txtrect = 'Measure area';
text(0.44*imWidth,0.4*imHeight,txtrect,'Color','y','fontsize',12);

txtcamera = 'Real time Camera Image';
text(0.36*imWidth,0.05*imHeight,txtcamera,'Color','w','fontsize',14);

txtimage_contrast=num2str(0);
txtimage_FFT=num2str(0);

%% 2) Contrast measurement routine
% It is possible to repeat this part to update the contrast and fft area
% calculation results on the camera preview screen
%
% It is not 100% real-time measurements, but it works pretty fast which
% will make the setting the lens position easier

% Clear text on the camera preview (contrast and FFT results)
delete(txtimage_contrast);
delete(txtimage_FFT);

% Save a camera image
start(vid);
image = getdata(vid);
imagesize = size(image);
Xpixel = imagesize(1);
Ypixel = imagesize(2);

% Cut image of the DMD (a,b,c,d values are from above)
a = 0.46*Xpixel;
b = 0.54*Xpixel;
c = 0.45*Ypixel;
d = 0.55*Ypixel;
imagecrop = image(a:b,c:d);
[Ypixel_crop Xpixel_crop] = size(imagecrop);

% figure(2); imshow(imagecrop); title('Cropped Test Image');

% Draw 25% / 50% / 75% position of the cropped image
imagecrop_25 = imagecrop(round(0.25*Ypixel_crop),:);
imagecrop_50 = imagecrop(round(0.50*Ypixel_crop),:);
imagecrop_75 = imagecrop(round(0.75*Ypixel_crop),:);
imagecrop_avg = mean([imagecrop_25;imagecrop_50;imagecrop_75]);
% imagecrop_std = std([imagecrop_25;imagecrop_50;imagecrop_75]);



%  3) Contrast measure 
white = max(imagecrop_avg);
black = min(imagecrop_avg);
contrast = (white-black)/(white+black)

% figure(3); hold on;
% plot(1:Xpixel_crop,imagecrop_25(:),'k-');
% plot(1:Xpixel_crop,imagecrop_50(:),'g-');
% plot(1:Xpixel_crop,imagecrop_75(:),'b-');
% plot(1:Xpixel_crop,imagecrop_avg(:),'r--','LineWidth',2);
% title('Line spread function');
% xlabel('Pixels (Horizontal)');
% ylabel('Intensity');
% legend('25%','50%','75%','Avg');



% 4) FFT area measure
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
% figure(4); hold on;
% plot(Psi,FFT_25(:),'k-');
% plot(Psi,FFT_50(:),'g-');
% plot(Psi,FFT_75(:),'b-');
% plot(Psi,FFT_avg(:),'r--');
% xlabel('Frequency distribution');
% ylabel('');
% legend('25%','50%','75%','Avg');

sumFFT_25 = sum(FFT_25);
sumFFT_50 = sum(FFT_50);
sumFFT_75 = sum(FFT_75);
sumFFT_avg = sum(FFT_avg)

% Close the cropped DMD image
% close(figure(2));

% Display measured contrast and FFT area results on the camera preview
txtcamera_contrast = append('Contrast:  ',num2str(contrast));
txtimage_contrast=text(1.3*imWidth*markerindex,0.7*imHeight*markerindex,txtcamera_contrast,'Color','w');

txtcamera_FFT = append('FFT area:  ',num2str(sumFFT_avg));
txtimage_FFT = text(1.3*imWidth*markerindex,0.8*imHeight*markerindex,txtcamera_FFT,'Color','w');
