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
clear; close all;

%% Open camera.
rectRatioWidth = 0.1;
rectRatioHeight = 0.08;
vid = OpenCamera('rectRatioHeight',rectRatioHeight,'rectRatioWidth',rectRatioWidth);

%% Capture image and calculate contrast (SACCSFA).
%
% Set which channel and spatial frequency to test.
numChannel = 15;
cyclesPerDegOptions = [3 6 9 12 18];
numViewingMedia = 1;

% Make a loop to measure for all spatial frequencies.
nSFs = length(cyclesPerDegOptions);
for ss = 1:nSFs
    % Set the target spatial frequency.
    cyclesPerDeg = cyclesPerDegOptions(ss);
    
    % Get a key stroke before starting each meausrement.
    fprintf('Press any to start measurement - (Ch %d) and (%d cpd) \n',numChannel,cyclesPerDeg);
    pause;
    
    % Clear text on the camera preview. This would make sort of real-time
    % measurement by updating the numbers.
    fig = gcf;
    textObjects = findall(fig,'Type','text');
    delete(textObjects);
    
    % Set which viewing media to use either SACCSFA or printed target for
    % camera MTF.
    switch numViewingMedia
        case 1
            viewingMedia = 'SACCSFA';
        case 2
            viewingMedia = 'Print';
    end
    
    % Set it differently over the channel. This value was found when we turn on
    % a single channel of combi-LED in input setting of 0.5 intensity (0-1).
    switch viewingMedia
        case 'SACCSFA'
            exposureTimePerChannel = [215000 0 80000 0 160000 990000 180000 100000 0 170000 0 120000 115000 0 120000 0];
        case 'Print'
            exposureTimePerChannel = [27000 45000 40000 58000 21000 210000 47000 47000];
    end
    exposureTime = exposureTimePerChannel(numChannel);
    
    % Set min peak distance between peaks for calculating spatial frequency.
    switch cyclesPerDeg
        case 3
            minPeakDistancePixel = 40;
        case 6
            minPeakDistancePixel = 17;
        case 9
            minPeakDistancePixel = 5;
        case 12
            minPeakDistancePixel = 5;
        case 18
            minPeakDistancePixel = 5;
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
end


%% Open the combi-LED projector before to use.
%
% Open the projector
obj = CombiLEDcontrol();

%% Capture image and calculate contrast (PRINT).
%
% Set which channel and spatial frequency to test.
cyclesPerDeg = 3;
channelOptions = [1 2 3 4 5 6 7 8];
channelIntensity = 0.5;
maxChannelIntensity = 1;
numViewingMedia = 2;

% Make a loop to measure for all spatial frequencies.
nChannels = length(channelOptions);
for ss = 1:nChannels
    % Set the target channel.
    channelSettings = [0,0,0,0,0,0,0,0];
    numChannel = channelOptions(ss);
    
    % Turn on the target channel.
    channelSettings(numChannel) = channelIntensity;
    
    % Send the settings
    obj.setPrimaries(channelSettings);
    
    % Get a key stroke before starting each meausrement.
    fprintf('Press any to start measurement - (Ch %d: %.1f/%.1f) and (%d cpd) \n',numChannel,channelIntensity,maxChannelIntensity,cyclesPerDeg);
    pause;
    
    % Clear text on the camera preview. This would make sort of real-time
    % measurement by updating the numbers.
    fig = gcf;
    textObjects = findall(fig,'Type','text');
    delete(textObjects);
    
    % Set which viewing media to use either SACCSFA or printed target for
    % camera MTF.
    switch numViewingMedia
        case 1
            viewingMedia = 'SACCSFA';
        case 2
            viewingMedia = 'Print';
    end
    
    % Set it differently over the channel. This value was found when we turn on
    % a single channel of combi-LED in input setting of 0.5 intensity (0-1).
    switch viewingMedia
        case 'SACCSFA'
            exposureTimePerChannel = [215000 0 80000 0 160000 990000 180000 100000 0 170000 0 120000 115000 0 120000 0];
        case 'Print'
            exposureTimePerChannel = [27000 45000 40000 58000 21000 210000 47000 47000];
    end
    exposureTime = exposureTimePerChannel(numChannel);
    
    % Set min peak distance between peaks for calculating spatial frequency.
    switch cyclesPerDeg
        case 3
            minPeakDistancePixel = 40;
        case 6
            minPeakDistancePixel = 17;
        case 9
            minPeakDistancePixel = 5;
        case 12
            minPeakDistancePixel = 5;
        case 18
            minPeakDistancePixel = 5;
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
end

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
