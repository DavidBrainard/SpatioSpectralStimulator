% SetCameraCenterPoint
%
% This shows a camera preview realtime to place the camera in parallel to
% the cage rod system for SACC project. To utilize this code, the USB
% camera should be connected. It basically uses the Image acquisition
% toolbox, and the add-ons 'Gentl' and 'Gig' should be also installed to
% identify the camera.

% History:
%    11/30/21 smo    Cleaned up and tried to write it in the format.
%    12/03/21 smo    Removed the finding centered point part and make it as
%                    a separate code (see FindCircleImage.m)

%% Initialzie.
clear; close all;

%% Identify the camera and setups.
%
% Load the camera information.
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview.
% idxWindowSize decides the size of the camera realtime preview [0-1].
idxWindowSize = 0.25;
vidResResize = vidRes * idxWindowSize;

hFig = figure('Units', 'pixels', 'Position', [100 100 vidResResize(1) vidResResize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidResResize(1) vidResResize(2)]);
hImage = image(zeros(vidRes(2), vidRes(1), nBands));

% Place the cross marker at the center of the preview screen.
imWidth = vidRes(1);
imHeight = vidRes(2);
numBands = vid.NumberOfBands;
markerIndex = 0.5;

% Display real-time camera preview screen with the center point marked.
preview(vid, hImage);
centerPoint = line(hAxes, round([markerIndex*imWidth, markerIndex*imWidth]),round([markerIndex*imHeight, markerIndex*imHeight]),'Marker','x','MarkerSize',15,'color','r','LineWidth',1,'LineStyle','none');
