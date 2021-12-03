% FindCircleImage
%
% This searches circles on a test image by desired sizes and positions.
% It is usefull to align the cage rod system and the camera for SACC project.

% History:
%    12/03/21 smo   Pulled out this part and made it as a separate code.
%                   (See also SetCameraCenterPoint.m)


%% Initialize.
clear all; close all;

%% Capture a new image if you want.
% 
% You can use this code either newly measuring the image or simply read the
% existing image to analyze.
%
% This part captures and saves the image and finds the circle on the image
% to decide it as the centered point. 

MEASURE = false;
if (MEASURE)
    % Set up the camera settings. Usually, you don't need to change this
    % part.
    %
    % FramesPerTrigger can be set within 1-255
    vid = videoinput('gentl', 1, 'Mono8');
    src = getselectedsource(vid);
    vid.FramesPerTrigger = 1; 
    src.AutoExposureTimeLowerLimit = 48;
    src.AutoExposureTimeUpperLimit = 10000;
    src.AutoFunctionROIWidth = 2780;
    
    % Start the video and save the captured image to analyze.
    start(vid);
    vidData = getdata(vid);
    testFileName = 'Test';
    imwrite(vidData,append(testFileName,'.tiff'));
end

%% Load an image and find circles.
if (~MEASURE)
    testFileName = 'Screen';
end
imageOriginal = imread(append(testFileName,'.tiff'),'tif');
[Ypixel Xpixel] = size(imageOriginal);

% Set an original image.
figure; clf;
imshow(imageOriginal); 
title('Original Image');   

% Find size of circles we want. First, put a randome line on the image
% and measure the radius of the target circle image.
rulerOnImage = imdistline; 
imageFiltered = imfilter(imageOriginal,ones(5)/25);

% Set searching parameters. Usually change the min and max size of the
% circles.
% 
% Note that the function 'imfindcircles' with small raidus value (less than
% or equal to 5, so diameter 10) would be limited in its accuracy.
minSizeCircleDiameter = 10;
maxSizeCircleDiameter = 40;
searchSensitivity = 0.95;
searchEdgeThreshold = 0.04;

% Find circles and mark on the image.
[centers, radii] = imfindcircles(imageFiltered, [minSizeCircleDiameter maxSizeCircleDiameter]/2, ...
                                'ObjectPolarity','dark','Sensitivity',searchSensitivity,'EdgeThreshold',searchEdgeThreshold);
h = viscircles(centers,radii);
    
% Find the circles with desired sizes at certain positions. We will set the
% position very narrow to sort out the circles at the center of the image.
figure; clf; hold on;
imshow(imageFiltered); 
title('Image with targeted points'); 

% Limit the range for inclduing the circles on the image.
minRange_XPixel = 0.495;
maxRange_XPixel = 0.505;
minRange_YPixel = 0.495;
maxRange_YPixel = 0.505;
    
idxTarget_XPixel = find(centers(:,1) > Xpixel*minRange_XPixel & centers(:,1) < Xpixel*maxRange_XPixel); 
idxTarget_YPixel = find(centers(:,2) > Ypixel*minRange_YPixel & centers(:,2) < Ypixel*maxRange_YPixel);

idxTarget = intersect(idxTarget_XPixel,idxTarget_YPixel);  
radiiTarget = radii(idxTarget);
centersTarget = centers(idxTarget,:);
hTarget = viscircles(centersTarget,radiiTarget);
    
% Find the coordinates of the centered point here.
centeredPosition = [centersTarget(1)/Xpixel centersTarget(2)/Ypixel];

% Add centered point coordinates on the image.
addOnTextImage = num2str(round(centeredPosition,4));
text(centersTarget(1),centersTarget(2),addOnTextImage,'HorizontalAlignment','left','color','r');
