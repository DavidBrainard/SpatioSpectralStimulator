% CameraContrastRealtime.
%
% This captures the real time camera image and calculate the contrast from
% the targeted area of the image.

% History:
%    06/13/23    smo     - Cleaned up from the old script.
%    08/31/23    smo     - Moved to project repository as Trombone laptop
%                          is now able to use git.

%% Open camera preview.
%
% Initialize.
clear; close all; clc;

% Load the camera info.
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview
% The screen ratio to original. 0.25 means 25% raw image (type 0 to 1).
windowsize = 0.25;
vidRes_resize = vidRes.*windowsize;
hFig = figure('Units', 'pixels', 'Position', [100 100 vidRes_resize(1) vidRes_resize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidRes_resize(1) vidRes_resize(2)]);
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

% Info about the marker position at the center.
imWidth = vidRes(1);
imHeight = vidRes(2);
numBands = vid.NumberOfBands;
markerindex = 0.5;

% Show camera preview here.
preview(vid, hImage)
hLine = line(hAxes, round([markerindex*imWidth, markerindex*imWidth]),...
    round([markerindex*imHeight, markerindex*imHeight]),...
    'Marker','+','MarkerSize',30,'color','r','LineWidth',1,'LineStyle','none');

% Set the targeted area of the image.
% a = rect starting x-coordinate, b =rect starting y-coordinate,
% c = rect width, d = rect height.
ratioWidth = 0.1;
ratioHeight = 0.08;
a = (0.5-ratioWidth/2)*imWidth;
b = (0.5-ratioHeight/2)*imHeight;
c = ratioWidth * imWidth;
d = ratioHeight * imHeight;

% Show the target area in a rectangle.
rectangle('Position',[a,b,c,d],'Curvature',[0,0],'LineWidth',1,'LineStyle','--','edgecolor','y')
txtrect = 'Measure area';
text(0.44*imWidth,0.43*imHeight,txtrect,'Color','y','fontsize',12);

% Set the title of the image.
txtcamera = 'Real time Camera Image';
text(0.36*imWidth,0.05*imHeight,txtcamera,'Color','w','fontsize',14);

%% Measure contrast.
%
% Repeat this part to update the contrast calculation results on the camera
% preview screen. It is not 100% real-time measurements, but works pretty
% fast.
%
% Clear text on the camera preview. This would make sort of real-time
% measurement by updating the numbers.
fig = gcf;
textObjects = findall(fig,'Type','text');
delete(textObjects);

% Save a camera image. We will calculate the contrast from still image.
start(vid);
image = getdata(vid);
imagesize = size(image);
pixelHeight = imagesize(1);
pixelWidth = imagesize(2);

% Save the raw image here for saving it later.
imageRaw = image(1:pixelHeight,1:pixelWidth);

% Get the targeted area of the image. 
a = (0.5-ratioHeight/2)*pixelHeight;
b = (0.5+ratioHeight/2)*pixelHeight;
c = (0.5-ratioWidth/2)*pixelWidth;
d = (0.5+ratioWidth/2)*pixelWidth;
imageCrop = image(a:b,c:d);
[Ypixel_crop Xpixel_crop] = size(imageCrop);

% Crop the targeted area in the image.
imagecrop_25 = imageCrop(round(0.25*Ypixel_crop),:);
imagecrop_50 = imageCrop(round(0.50*Ypixel_crop),:);
imagecrop_75 = imageCrop(round(0.75*Ypixel_crop),:);
imagecrop_avg = mean([imagecrop_25;imagecrop_50;imagecrop_75]);

% Calculate contrast.
white = max(imagecrop_avg);
black = min(imagecrop_avg);
contrast = (white-black)/(white+black);
fprintf('Contrast = (%.2f) \n', contrast);

% Show contrast on the camera preview.
textContrast = append('Contrast:  ',num2str(round(contrast,2)));
text(1.05*imWidth*markerindex, 1*imHeight*markerindex, textContrast, 'Color','w');

% Show spatial frequency.
%
% Set minimum peak distance here for not having multiple peaks at one peak
% point. Followings are recommendation over different spatial frequency.
% 
% minPeakDistance = 5 (9, 12 ,18 cpd), 17 (6 cpd), 40 (3 cpd).
minPeakDistance = 40;

% Show the peaks found. Visually check.
figPeak = figure; findpeaks(double(imagecrop_50),'MinPeakDistance',minPeakDistance)

% Calculate the spatial frequency here.
[~,peakIndex] = findpeaks(double(imagecrop_50),'MinPeakDistance',minPeakDistance);
numCycles = length(peakIndex);

% Get the horizontal size of cropped image in degrees.
pixelToInchHorizontal = 0.0367;
pixelToInchVertical = 0.0362;
physicalDistnaceRefInch = 370;
cropImageHalfInch = pixelToInchHorizontal * (Xpixel_crop/2);
sizeDegHorizontal = 2*(rad2deg(atan(cropImageHalfInch/physicalDistnaceRefInch)));

% When the cropped image was 309 x 166 pixels.
% sizeDegHorizontal = 1.7570;

% Calculate the cpd here.
cyclesPerDeg = numCycles/sizeDegHorizontal;
fprintf('Spatial frequency = (%.1f) \n', cyclesPerDeg);

% Add spatial frequency on the camera preview.
textCyclesPerDeg = append('Spatial frequency:  ',num2str(round(cyclesPerDeg,0)));
text(1.05*imWidth*markerindex, 1.1*imHeight*markerindex, textCyclesPerDeg,'Color','w');

% Save the image if you want.
SAVEIMAGE = false;
if (SAVEIMAGE)
    cyclesPerDegStr = 'single_';
    fileNameRawImage = append(cyclesPerDegStr,'raw_');
    fileNameCropImage = append(cyclesPerDegStr,'crop_');
    dayTimeStr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    imwrite(imageCrop,append(fileNameCropImage,dayTimeStr,'.tiff'));
    imwrite(imageRaw,append(fileNameRawImage,dayTimeStr,'.tiff'));
    disp('Image has been saved successfully!');
end
