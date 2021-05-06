%% Get an image using Basler camera 

%% Initialize
clear; close all;
imaqhwinfo % make sure 'gentl' and 'gige' are installed (available thru add-ons) % Type 'imaqhelp' for further info

%% Configure the camera in its full resolution
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1; % set within 1-255
src.AutoExposureTimeLowerLimit = 48;
src.AutoExposureTimeUpperLimit = 10000;
src.AutoFunctionROIWidth = 2780;
preview(vid); % real-time camera scene % To stop: stoppreview(vid);

%% Set a smaller preview window
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


%% Save the captured image as a desired image file format
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1; % set within 1-255
src.AutoExposureTimeLowerLimit = 48;
src.AutoExposureTimeUpperLimit = 10000;
src.AutoFunctionROIWidth = 2780;

start(vid);
viddata = getdata(vid);
imwrite(viddata,'Test_0506.tiff');

%% Preview Webcam Live Video Stream with Custom Window
% 
% fig = figure('NumberTitle','off','MenuBar','none');
% fig.Name = 'My Camera';
% ax = axes(fig); 
% start(vid); % take a single shot to fit the frame
% frame = ans;
% im = image(ax,zeros(size(frame),'uint8')); 
% axis(ax,'image');
% 
% preview(vid,im)

%% Others

% % Initialize
% clear; close all;
% 
% webcamlist % show all connected cam list
% cam = webcam; % use the first camera returned by 'webcamlist'
% 
% preview(cam);
% 
% cam.AvailableResolutions % Available resolution settings
% cam.Resolution = '1280x720';
% 
% img = snapshot(cam); % Take a snap
% imshow(img);
% 
% % Loop 30 times
% for frames = 1:30
%     img = snapshot(cam);
%     imsho(img)
% end
