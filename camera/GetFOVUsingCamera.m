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
%    04/27/22   smo   Added the complete data of the values for both
%                     vertical and horizontal references.
%    05/10/22   smo   Added an option to calculate the FOV by reading the
%                     scale on the referece (ruler) directly.

%% Initialize.
clear; close all;

%% Set parameters.
VERBOSE = true;
CALCULATIONMETHOD = 'direct';

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
% Reference size in inch. We meausred both horizontal and vertical
% measurements.
%
% The length is the end-to-end of either vertical or horizontal camera
% captured image.
DATE = '0427';

switch DATE
    case '0426'
        % 04/26 data.
        refHorizontalInch = 41.75;
        refVerticalInch = 85.50;
    case '0427'
        % 04/27 data.
        refHorizontalInch = 113.40;
        refVerticalInch = 74.80;
end

% Reference size in pixel. This is the same with the camera resolution.
refVerticalPixel = imgSize(1);
refHorizontalPixel = imgSize(2);

% Get inch per pixel. You can choose either vertical or horizontal sides to
% use.
whichSideRef = 'horizontal';

switch whichSideRef
    case 'vertical'
        inchPerPixel = refVerticalInch/refVerticalPixel;
        % Distance between the camera and the reference ruler in inch.
        distanceInch = 429;
        if (DATE == '0427')
            distanceInch = 370;
        end
    case 'horizontal'
        inchPerPixel = refHorizontalInch/refHorizontalPixel;
        distanceInch = 140;
        if (DATE == '0427')
            distanceInch = 370;
        end
end

% Target image length in inch.
if (strcmp(CALCULATIONMETHOD,'indirect'))
    targetImgHorizontalLeftInch = targetImgHorizontalLeftPixel * inchPerPixel;
    targetImgHorizontalRightInch = targetImgHorizontalRightPixel * inchPerPixel;
    
    targetImgVerticalTopInch = targetImgVerticalLeftPixel * inchPerPixel;
    targetImgVerticalBottomInch = targetImgVerticalRightPixel * inchPerPixel;
end

%% Read the scale directly from the reference image.
%
% Read the reference image here.
refDir = 'horizontal';
refDate = '2022-04-27';
fileType = '.tiff';
fileName = append(refDir,'_',refDate,fileType);
refImg = imread(fileName);

% Plot it.
figure;
imshow(refImg); hold on;

% Mark the DMD outline on the reference image.
rectangle('Position', targetImgCords, 'EdgeColor','r', 'LineWidth', 3);

% Mark the center of image.
plot(imgCenter(2),imgCenter(1),'yo','MarkerFaceColor','y','MarkerSize',12);
line([0 imgSize(2)],[imgCenter(1) imgCenter(1)],'color','yellow','linewidth',2)
line([imgCenter(2) imgCenter(2)],[0 imgSize(1)],'color','yellow','linewidth',2)
legend('Center of the image','fontsize',15);

% Target image length in inch. You can read the scale on the ruler directly
% in the reference image.
if (strcmp(CALCULATIONMETHOD,'direct'))
    targetImgHorizontalLeftInch = 74 - 32.3;
    targetImgHorizontalRightInch = 121.4 - 74;
    
    targetImgVerticalTopInch = 46.5 - 18.6;
    targetImgVerticalBottomInch = 68.6 - 46.5;
end

%% Calculate the FOV of the target in degrees.
%
% Horizontal.
targetImgHorizontalDegLeft = rad2deg(atan(targetImgHorizontalLeftInch/distanceInch));
targetImgHorizontalDegRight = rad2deg(atan(targetImgHorizontalRightInch/distanceInch));
targetImgHorizontalDeg = targetImgHorizontalDegLeft + targetImgHorizontalDegRight;

% Vertical.
targetImgVerticalDegLeft = rad2deg(atan(targetImgVerticalTopInch/distanceInch));
targetImgVerticalDegRight = rad2deg(atan(targetImgVerticalBottomInch/distanceInch));
targetImgVerticalDeg = targetImgVerticalDegLeft + targetImgVerticalDegRight;

% Combined.
targetImgDeg = [targetImgHorizontalDeg targetImgVerticalDeg]
