function [vid] = OpenCamera(options)
% Open the camera and show the preview.
%
% Syntax:
%    [vid] = OpenCamera()
%
% Description:
%    This is to open and make it ready to use the camera to capture the
%    image. The CCD camera should be connected to the computer in advance.
%    It will show the camera preview which shows the real-time scene that
%    camera captures.
%
%    On the image of the preivew, we will show the marker at the center and
%    also rectangle where we eventually use to measure contrast. The size
%    of the rectangle can be adjusted.
%
% Inputs:
%    N/A
%
% Outputs:
%    vid                      - Video object to control.
%
% Optional key/value pairs:
%    windowSize               - Window size of the preview to show. Default
%                               to 0.25. It's proportion to the full
%                               resolution, so 0.25 means its quarter size
%                               of the full resolution.
%    rectRatioWidth           - Define the width of the rectangle in
%                               proportion to the full resolution. Default
%                               to 0.1.
%    rectRatioHeight          - Define the height of the rectangle in
%                               proportion to the full resolution. Default
%                               to 0.08.
%
% See also:
%    CameraContrastRealtime

% History:
%    11/29/23       smo       - Wrote it

%% Set variables.
arguments
    options.windowSize (1,1) = 0.25
    options.rectRatioWidth (1,1) = 0.1
    options.rectRatioHeight (1,1) = 0.08 
end

%% Open camera preview.
%
% Load the camera info.
vid = videoinput('gentl', 1, 'Mono8');
vidRes = get(vid, 'VideoResolution');
nBands = get(vid, 'NumberOfBands');

% Set the preferred window size to preview
% The screen ratio to original. 0.25 means 25% raw image (type 0 to 1).
vidRes_resize = vidRes.* options.windowSize;
hFig = figure('Units', 'pixels', 'Position', [100 100 vidRes_resize(1) vidRes_resize(2)]);
hAxes = axes('Units', 'pixels', 'Position', [10 10 vidRes_resize(1) vidRes_resize(2)]);
hImage = image(zeros(vidRes(2), vidRes(1), nBands));

% Show the camera preview here.
preview(vid, hImage)

% Show the marker at the center of the preview.
imWidth = vidRes(1);
imHeight = vidRes(2);
markerindex = 0.5;
hLine = line(hAxes, round([markerindex*imWidth, markerindex*imWidth]),...
    round([markerindex*imHeight, markerindex*imHeight]),...
    'Marker','+','MarkerSize',30,'color','r','LineWidth',1,'LineStyle','none');

% Set the size of the rectangle to draw on the preview image.
%
% a = rect starting x-coordinate
% b = rect starting y-coordinate
% c = rect width
% d = rect height
a = (0.5 - options.rectRatioWidth/2) * imWidth;
b = (0.5 - options.rectRatioHeight/2) * imHeight;
c = options.rectRatioWidth * imWidth;
d = options.rectRatioHeight * imHeight;

% Draw the rectangle here. It will be shown in the yellow lines.
rectangle('Position',[a,b,c,d],'Curvature',[0,0],'LineWidth',1,'LineStyle','--','edgecolor','y')

% Add text on the preview.
textRect = 'Measure area';
text(0.44*imWidth,0.43*imHeight,textRect,'Color','y','fontsize',12);

% Show the title of the image.
txtcamera = 'Real time Camera Image';
text(0.36*imWidth,0.05*imHeight,txtcamera,'Color','w','fontsize',14);

% Print out if all worked.
fprintf('Camera is connected and working now! \n Rect area = %.0f x %.0f (pixels) \n', c, d);

end
