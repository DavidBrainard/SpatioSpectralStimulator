% CameraContrastRealTime
%
% This is written to place the Lens 1 on the cage rod system for SACC
% project. What it does is to show the real-time camera view and calculates
% the contrast of the image. Test image should be vertical stripe image (we
% used 20 pixels/bar image).

% History:
%    12/06/21  smo   Started on cleaning.

%% Initialize.
clear; close all; 

%% Show the real-time camera image.
% 
% Load the camera info
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview
windowSize = 0.25; 
vidResResize = vidRes .* windowSize;

hFig = figure('Units', 'pixels', 'Position', [100 100 vidResResize(1) vidResResize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidResResize(1) vidResResize(2)]);
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

% Info about the marker position.
% markerIndex sets the centered point.
imWidth = vidRes(1);
imHeight = vidRes(2);
numBands = vid.NumberOfBands;
markerIndex = 0.5; 

% Camera preview with the centered point marked
figure; clf;
preview(vid, hImage);

% Set the marker here.
hLine = line(hAxes, round([markerIndex*imWidth, markerIndex*imWidth]),round([markerIndex*imHeight, markerIndex*imHeight]),...
        'Marker','+','MarkerSize',30,'color','r','LineWidth',1,'LineStyle','none');

% Set the rectangle here.
fromX = 0.45 * imWidth; 
fromY = 0.46 * imHeight; 
rectWidth = 309; 
rectHeight = 166; 
rectangle('Position',[fromX,fromY,rectWidth,rectHeight],'Curvature',[0,0],'LineWidth',1,'LineStyle','--','edgecolor','y')
txtRect = 'Measure area';
text(0.44*imWidth, 0.4*imHeight, txtRect, 'Color', 'y', 'fontsize', 12);

txtCamera = 'Real time Camera Image';
text(0.36*imWidth, 0.05*imHeight, txtCamera, 'Color', 'w', 'fontsize', 14);

% Set the initial value to zero.
txtContrast = num2str(0);
txtFFT = num2str(0);

%% Contrast measurement routine.
%
% It is possible to repeat this part to update the contrast and FFT area
% calculation results on the camera preview screen.
%
% It is not 100% real-time measurements, but it works pretty fast which
% will make placing the lens 1 position easier.

% Clear text on the camera preview. It deletes the previous calculation
% results.
delete(txtContrast);
delete(txtFFT);

% Get the current camera image.
start(vid);
image = getdata(vid);
imageSize = size(image);
screenXpixel = imageSize(1);
screenYpixel = imageSize(2);

% Clip the part of the screen image to calculate the contrast.
fromX = 0.46 * screenXpixel;
fromY = 0.54 * screenXpixel;
rectWidth = 0.45 * screenYpixel;
rectHeight = 0.55 * screenYpixel;
imageCrop = image(fromX:fromY,rectWidth:rectHeight);
[YpixelCrop XpixelCrop] = size(imageCrop);

% We will use the average of the 25% / 50% / 75% positions of the cropped image. 
imageCrop25 = imageCrop(round(0.25*YpixelCrop),:);
imageCrop50 = imageCrop(round(0.50*YpixelCrop),:);
imageCrop75 = imageCrop(round(0.75*YpixelCrop),:);
imageCropAvg = mean([imageCrop25;imageCrop50;imageCrop75]);

% Calculate contrast here.
whiteCropImage = max(imageCropAvg);
blackCropImage = min(imageCropAvg);
contrastCropImage = (whiteCropImage-blackCropImage) / (whiteCropImage+blackCropImage);

% Calculate FFT here. This is just additional measure to check the
% contrast.
FFT25  = abs(fftshift(fft(imageCrop25)));
FFT50  = abs(fftshift(fft(imageCrop50)));
FFT75  = abs(fftshift(fft(imageCrop75)));
FFTAvg = abs(fftshift(fft(imageCropAvg)));

% Normalize.
FFT25  = FFT25/max(FFT25);
FFT50  = FFT50/max(FFT50);
FFT75  = FFT75/max(FFT75);
FFTAvg = FFTAvg/max(FFTAvg);

% Set x-axis of the FFT results (frequecy).
% This part needs to be updated
Size = length(FFT25);
spacing = 20;% spacing between data points 
fsx = abs(1/spacing); % turn into sampling frequency
fromX = linspace(-Size/2,Size/2,Size); % form scale for conversion based on frequency bins
conversionx = fsx./Size; % conversion factor for frequency bin units to frequency (unit^-1)
Psi = fromX.*conversionx; % frequency (unit^-1)

sumFFT25  = sum(FFT25);
sumFFT50  = sum(FFT50);
sumFFT75  = sum(FFT75);
sumFFTAvg = sum(FFTAvg);

% Display calculation results on the camera preview.
txtContrast = append('Contrast: ', num2str(contrastCropImage));
txtContrast = text(1.3*imWidth*markerIndex, 0.7*imHeight*markerIndex, txtContrast, 'Color', 'w');

txtFFT = append('FFT area: ', num2str(sumFFTAvg));
txtFFT = text(1.3*imWidth*markerIndex, 0.8*imHeight*markerIndex, txtFFT, 'Color', 'w');
