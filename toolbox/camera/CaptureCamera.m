function [image] = CaptureCamera(vid,options)
% Capture an image using a camera.
%
% Syntax:
%    [image] = CaptureCamera()
%
% Description:
%    This captures an image of the scene using a camera connected to the
%    computer.
%
% Inputs:
%    N/A
%
% Outputs:
%    vid                      - Video object to control.
%
% Optional key/value pairs:
%    dd
%
% See also:
%    OpenCamera, CameraContrastRealtime

% History:
%    11/30/23       smo       - Wrote it

%% Set variables.
arguments
    vid
    options.exposureTime (1,1) = 10000
    options.rectRatioHeight (1,1) = 0.08
    options.rectRatioWidth (1,1) = 0.1
    options.minPeakDistancePixel (1,1) = 10
    options.saveImageDir = []
    options.saveImageFilename
end

%% Control exposure time.
src = getselectedsource(vid);
src.ExposureTime = options.exposureTime;

%% Save a camera image. We will calculate the contrast from still image.
start(vid);
image = getdata(vid);
imageSize = size(image);
imagePixelHeight = imageSize(1);
imagePixelWidth = imageSize(2);

% Save the raw image here for saving it later.
rawImage = image(1:imagePixelHeight, 1:imagePixelWidth);

% Get the size of the targeted area within the image.
a = (0.5-options.rectRatioHeight/2)*imagePixelHeight;
b = (0.5+options.rectRatioHeight/2)*imagePixelHeight;
c = (0.5-options.rectRatioWidth/2)*imagePixelWidth;
d = (0.5+options.rectRatioWidth/2)*imagePixelWidth;

% Get the image here.
cropImage = image(a:b,c:d);
[Ypixel_crop Xpixel_crop] = size(cropImage);

% We calculate the mean of three (25%, 50%, 75%) intensity profiles.
imagecrop_25 = cropImage(round(0.25*Ypixel_crop),:);
imagecrop_50 = cropImage(round(0.50*Ypixel_crop),:);
imagecrop_75 = cropImage(round(0.75*Ypixel_crop),:);
imagecrop_avg = mean([imagecrop_25;imagecrop_50;imagecrop_75]);

% Calculate contrast from the averaged intensity profile. The contrast here
% is simply picking the highest and the lowest for white and black,
% responsively. This is just for quick checking. For SACC project, we
% calculate contrasts properly for analysis.
white = max(imagecrop_avg);
black = min(imagecrop_avg);
contrast = (white-black)/(white+black);
fprintf('Contrast = (%.2f) \n', contrast);

% Show contrast on the camera preview. This would make real-time ish
% measurement.
markerIndex = 0.5;
vidRes = get(vid, 'VideoResolution');
imWidth = vidRes(1);
imHeight = vidRes(2);
textContrast = sprintf('Contrast: %.2f',contrast);
text(1.05*imWidth*markerIndex, 1*imHeight*markerIndex, textContrast, 'Color', 'w');

%% Calculate spatial frequency of the image.
%
% Find the peaks in the intensity profile. We recommend to visually check
% the number of peaks found correctly, then proceed to calculate spatial
% frequency.
figure;
findpeaks(double(imagecrop_50),'MinPeakDistance',options.minPeakDistancePixel)

% Calculate the spatial frequency here.
[~,peakIndex] = findpeaks(double(imagecrop_50),'MinPeakDistance',options.minPeakDistancePixel);
numCycles = length(peakIndex);

% Get the horizontal size of cropped image in degrees. These values were
% found from the measurement, which would not be changed for SACC project.
% We will calculate the horizontal size in degree to compute spatial
% frequency.
%
% We first calculate the degree of the half angle, then multiply by 2 to
% get the visual angle of the scene.
pixelToInchHorizontal = 0.0367;
physicalDistnaceRefInch = 370;
cropImageHalfInch = pixelToInchHorizontal * (Xpixel_crop/2);
sizeDegHorizontal = 2*(rad2deg(atan(cropImageHalfInch/physicalDistnaceRefInch)));

% Calculate spatial frequency here.
cyclesPerDeg = numCycles/sizeDegHorizontal;
fprintf('Spatial frequency = (%.1f) \n', cyclesPerDeg);

% Add spatial frequency on the camera preview.
textCyclesPerDeg = sprintf('Spatial frequency: %.1f',cyclesPerDeg);
text(1.05*imWidth*markerIndex, 1.1*imHeight*markerIndex, textCyclesPerDeg, 'Color', 'w');

% Save the image if you want.
if ~isempty(options.saveImageDir)
    testfileDir = options.saveImageDir;
    
    % Make a new directory if there is no directory with the name.
    if ~isdir(testfileDir)
        mkdir(testfileDir);
    end
    
    % Get the date for the file name.
    dayTimeStr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    
    % Save the raw image.
    fileNameRawImage = fullfile(testfileDir,append(options.saveImageFilename,'_raw_'));
    imwrite(rawImage,append(fileNameRawImage,dayTimeStr,'.tiff'));
    
    % Save the cropped image.
    fileNameCropImage = fullfile(testfileDir,append(options.saveImageFilename,'_crop_'));
    imwrite(cropImage,append(fileNameCropImage,dayTimeStr,'.tiff'));
    
    disp('Image has been saved successfully!');
end

end
