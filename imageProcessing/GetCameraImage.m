%% Get an image using Basler camera 

% *** Before to start ***
% 1) Make sure the MATLAB is started with sudo command (cf. sudo matlab on terminal)


%% Initialize% 1
clear; close all; clc;
imaqhwinfo % make sure 'gentl' and 'gige' are installed (available thru add-ons) % Type 'imaqhelp' for further info

%% Configure the camera in its full resolution
% vid = videoinput('gentl', 1, 'Mono8');
% src =imaqhwinfo.getselectedsourvid = videoinput('gentl', 1, 'Mono8');
% vid.FramesPerTrigger = 1; % set within 1-255
% src.AutoExposureTimeLowerLimit = 48;
% src.AutoExposureTimeUpperLimit = 10000;
% src.AutoFunctionROIWidth = 2780;
% preview(vid); % real-time camera scene % To stop: stoppreview(vid);

%% Set a smaller preview window
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

windowsizeindex = 0.25;
vidRes_resize = vidRes.*windowsizeindex;

figure('Units', 'pixels', 'Position', [100 100 vidRes_resize(1) vidRes_resize(2)]);
axes('Units', 'pixels', 'Position', [10 10 vidRes_resize(1) vidRes_resize(2)]);
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

preview(vid, hImage);

%% Save the captured image as a desired image file format
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1; % set within 1-255
src.AutoExposureTimeLowerLimit = 48;
src.AutoExposureTimeUpperLimit = 10000;
src.AutoFunctionROIWidth = 2780;

start(vid);
viddata = getdata(vid);
imwrite(viddata,'LCPA1_1.tiff');

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

%% Others using 'Webcam support toolbox'

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
% Loop 30 times
% for frames = 1:30
%     img = snapshot(cam);
%     imsho(img)
% end
