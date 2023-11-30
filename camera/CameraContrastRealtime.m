% CameraContrastRealtime.
%
% This captures the real time camera image and calculate the contrast from
% the targeted area of the image.

% History:
%    06/13/23    smo     - Cleaned up from the old script.
%    08/31/23    smo     - Moved to project repository as Trombone laptop
%                          is now able to use git.
%    11/29/23    smo     - Parts were subsituted with function.

%% Initialize.
clear; close all; clc;

%% Open camera.
rectRatioWidth = 0.1;
rectRatioHeight = 0.08;
vid = OpenCamera('rectRatioHeight',rectRatioHeight,'rectRatioWidth',rectRatioWidth);

%% Capture image and calculate contrast.
%
% Clear text on the camera preview. This would make sort of real-time
% measurement by updating the numbers.
fig = gcf;
textObjects = findall(fig,'Type','text');
delete(textObjects);

% Set which channel and spatial frequency to test.
numChannel = 1;
cyclesPerDeg = 3;

% Set it differently over the channel. This value was found when we turn on
% a single channel of combi-LED in input setting of 0.5 intensity (0-1).
exposureTimePerChannel = [27000 45000 40000 58000 21000 210000 47000 47000];
exposureTime = exposureTimePerChannel(numChannel);

% Set which viewing media to use either SACCSFA or printed target for
% camera MTF.
numViewingMedia = 1;
switch numViewingMedia
    case 1
        viewingMedia = 'SACCSFA';
    case 2
        viewingMedia = 'Print';
end

% Get today's date.
date = datestr(now,'yyyy-mm-dd');

% We will add channel and spatial frequency in the file name.
channel = append('Ch',num2str(numChannel));
SF = append(num2str(cyclesPerDeg),'cpd');

% Set the directory to save the image.
saveImageDir = getpref('SpatioSpectralStimulator','SACCMaterials');
saveImageDir = fullfile(saveImageDir,'Camera','ChromaticAberration',viewingMedia,date,channel);

% Capture an image here and do analysis of it.
CaptureCamera(vid,'rectRatioHeight',rectRatioHeight,'rectRatioWidth',rectRatioHeight,...
    'exposureTime',exposureTime,'minPeakDistancePixel',minPeakDistancePixel,...
    'saveImagedir',saveImageDir,'saveImageFilename',SF);

%% Test for one image.
%
% Here we are not saving the image and just checking if everything works
% fine before running the routine.
CAPTUREONEIMAGE = false;

if (CAPTUREONEIMAGE)
    % Set variables.
    exposureTime = 10000;
    minPeakDistancePixel = 40;
    
    % Capture image here.
    CaptureCamera(vid,'rectRatioHeight',rectRatioHeight,'rectRatioWidth',rectRatioHeight,...
        'exposureTime',exposureTime,'minPeakDistancePixel',minPeakDistancePixel,...
        'saveImagedir',[],'saveImageFilename',[]);
end
