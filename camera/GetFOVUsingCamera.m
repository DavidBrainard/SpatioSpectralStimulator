% GetFOVUsingCamera.
%
% This calculates the FOV using camera meausured images. The aim of this
% was to calculate the FOV of the DMD of the projector for SACC project.
%
% The concept is to measure two different images, a reference image and a
% target DMD image, both with the infinity focus, then we can calculate the
% FOV of the DMD using the known distance/pixel information from the
% reference image.

% History
%    04/26/22   smo   Corrected the calculation method.

%% Initialize.
clear; close all;

%% Set parameters.
VERBOSE = true;

%% Load target image and get size in pixels.
imgFileName = 'infinity_DMD.tiff';
img = imread(imgFileName);
imgSize = size(img);
imgCenter = imgSize/2;

%% Find the target on the image.
%
% Make a binary image.
bwImg = imbinarize(img);

% Find both black and white regions.
% Note that BoundingBox contains [topleft X; topleft Y; X width; Y length].
imgStats = regionprops(bwImg);
targetImgCenter = imgStats.Centroid;
targetImgCords = imgStats.BoundingBox;

targetImgHorizontalLeftPixel  = abs(targetImgCords(1) - imgCenter(2));
targetImgHorizontalRightPixel = abs(targetImgCords(1)+targetImgCords(3) - imgCenter(2));

targetImgVerticalLeftPixel  = abs(targetImgCords(2) - imgCenter(1));
targetImgVerticalRightPixel = abs(targetImgCords(2)+targetImgCords(4) - imgCenter(1));

% Show the image and draw the detected rectangle on it.
if (VERBOSE)
    imshow(bwImg);
    hold on;
    for i = 1:numel(imgStats)
        rectangle('Position', imgStats(i).BoundingBox, ...
            'Linewidth', 3, 'EdgeColor', 'r', 'LineStyle', '--');
    end
    title('Target image with its shape found','fontsize',15);
end

%% Get reference image info.
%
% Reference size in inch.
%
% We meausred both horizontal and vertical measurements.
%
% For horizontal measurement, we used two rulers (24 inches each).
% For vertical measurement, we used elastic ruler (300 inches).
%
% The length typed here is the length for the end-to-end of either
% vertical or horizontal camera captured image.
refHorizontalInch = (23 + 6/16) + (18 + 6/16);
refVerticalInch = 94.5 - 10 + 1;

% Reference size in pixel. This is the same with the camera resolution.
refVerticalPixel = imgSize(1);
refHorizontalPixel = imgSize(2);

% Get inch per pixel. You can choose either vertical or horizontal sides to
% use.
whichSide = 'horizontal';

switch whichSide
    case 'vertical'
        inchPerPixel = refVerticalInch/refVerticalPixel;
        % Distance between the camera and the reference ruler in inch.
        distanceInch = 429;
    case'horizontal'
        inchPerPixel = refHorizontalInch/refHorizontalPixel;
        distanceInch = 140;
    otherwise
end

%% Calculate the FOV of the target in degrees.
%
% Horizontal.
targetImgHorizontalLeftInch = targetImgHorizontalLeftPixel * inchPerPixel;
targetImgHorizontalRightInch = targetImgHorizontalRightPixel * inchPerPixel;

targetImgHorizontalDegLeft = rad2deg(atan(targetImgHorizontalLeftInch/distanceInch));
targetImgHorizontalDegRight = rad2deg(atan(targetImgHorizontalRightInch/distanceInch));

targetImgHorizontalDeg = targetImgHorizontalDegLeft + targetImgHorizontalDegRight;

% Vertical.
targetImgVerticalLeftInch = targetImgVerticalLeftPixel * inchPerPixel;
targetImgVerticalRightInch = targetImgVerticalRightPixel * inchPerPixel;

targetImgVerticalDegLeft = rad2deg(atan(targetImgVerticalLeftInch/distanceInch));
targetImgVerticalDegRight = rad2deg(atan(targetImgVerticalRightInch/distanceInch));

targetImgVerticalDeg = targetImgVerticalDegLeft + targetImgVerticalDegRight;

% Combined.
targetImgDeg = [targetImgHorizontalDeg targetImgVerticalDeg]
