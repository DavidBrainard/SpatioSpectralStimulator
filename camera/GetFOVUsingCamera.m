% GetFOVUsingCamera.
%
% This calculates the FOV of the DMD.

%% Initialize.
clear; close all;

%% Get DMD size in pixels.
imgDMD = 'infinity_DMD.tiff';
img = imread(imgDMD);
imgSize = size(img);
imgCenter = imgSize/2;

maxImg = 255;
xPixelIndex = find(img(imgCenter(1),:) == maxImg);
xPixelLength = max(xPixelIndex) - min(xPixelIndex);

yPixelIndex = find(img(:,imgCenter(2)) == maxImg);
yPixelLength = max(yPixelIndex) - min(yPixelIndex);

%% Get ruler size in inch and deg.
%
% Distance between ruler and camera.
cameraDistanceRulerInch = 140;

% Ruler size in inch.
% Ruler is 24 inches long and the length adds up two different rulers.
horizontalLengthRulerInch = (23 + 6/16) + (18 + 6/16);
halfHorizontalLengthRulerInch = horizontalLengthRulerInch/2;

halfHorizontalRulerDeg = rad2deg(atan(halfHorizontalLengthRulerInch/cameraDistanceRulerInch));
horizontalRulerDeg = halfHorizontalRulerDeg * 2;

% Ruler size in pixel.
horizontalLengthRulerPixel = imgSize(2);

% Calculate deg per pixel.
degPerPixel = horizontalRulerDeg/horizontalLengthRulerPixel;

%% Get DMD size in deg.
dmdXPixelDeg = xPixelLength * degPerPixel;
dmdYPixelDeg = yPixelLength * degPerPixel;
