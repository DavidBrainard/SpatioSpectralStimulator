%% 1) Set a smaller preview window
clear; close all;

% Load the camera info
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview
windowsizeindex = 0.25;
vidRes_resize = vidRes.*windowsizeindex;

hFig = figure('Units', 'pixels', 'Position', [100 100 vidRes_resize(1) vidRes_resize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidRes_resize(1) vidRes_resize(2)]);
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

% Info about the marker position
imWidth = vidRes(1);
imHeight = vidRes(2);
numBands = vid.NumberOfBands;
markerindex = 0.5; % Set the centered point

% Camera preview with the centered point marked
preview(vid, hImage)
hLine = line(hAxes, round([markerindex*imWidth, markerindex*imWidth]),round([markerindex*imHeight, markerindex*imHeight]),'Marker','x','MarkerSize',15,'color','r','LineWidth',1,'LineStyle','none');

%% 2) Save the image and find the centered point
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1; % set within 1-255
src.AutoExposureTimeLowerLimit = 48;
src.AutoExposureTimeUpperLimit = 10000;
src.AutoFunctionROIWidth = 2780;

start(vid);
viddata = getdata(vid);
imwrite(viddata,'Test.tiff');

% Load an image and find targeted circle
image = imread('Test.tiff','tif');
[Ypixel Xpixel] = size(image);

figure(1); imshow(image); title('Original Image');

% Find size of circles we want
d = imdistline; % measure the radius of the target circle image

image_filtered = imfilter(image,ones(5)/25);

% LCPA1 image
[centers, radii] = imfindcircles(image_filtered,[10 20]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
% [10 18] for the position 1
% [14 22] for the position 2

h = viscircles(centers,radii);

% Find targeted circles 
figure(2); imshow(image_filtered); title('Image with targeted points'); hold on;

% LCPA1 centered point 
lowerlimit_X = 0.485;
upperlimit_X = 0.515;
lowerlimit_Y = 0.485;
upperlimit_Y = 0.515;

index_target_X = find(centers(:,1)>Xpixel*lowerlimit_X & centers(:,1)<Xpixel*upperlimit_X); % sort of handpicked targeted points...
index_target_Y = find(centers(:,2)>Ypixel*lowerlimit_Y & centers(:,2)<Ypixel*upperlimit_Y); 

index_target = intersect(index_target_X,index_target_Y);

radii_target = radii(index_target);
centers_target = centers(index_target,:); % [Xpixel Ypixel]
h_target = viscircles(centers_target,radii_target);

position_error = [centers_target(1)/Xpixel centers_target(2)/Ypixel]

txt = num2str(round(position_error,4));
text(centers_target(1),centers_target(2),txt,'HorizontalAlignment','left','color','r')
