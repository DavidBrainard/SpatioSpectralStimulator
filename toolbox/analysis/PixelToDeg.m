function [deg pixels_SACCSFA] = PixelToDeg(pixels_camera,options)
% Convert pixel unit to field of view (deg).
%
% Syntax:
%    [deg] = PixelToDeg(pixels)
%
% Description:
%    dd
%
% Inputs:
%    pixels                   - Target length in pixels to convert the unit.
%
% Outputs:
%    deg                      - Field of view corresponding the input pixel
%                               size in degrees.
%
% Optional key/value pairs:
%    verbose                  - Default to false. Controls plotting.
%
% See also:
%    GetFOVUsingCamera.

% History:
%    01/17/24  smo            - Wrote it.

%% Set variables.
arguments
    pixels_camera (1,1)
    options.dir = 'horizontal'
    options.verbose = 'false'
end

%% Get the reference image info.
%
% The reference image is the image taken far enough to focus with infinity
% camera focus setting. Here, the distance between the scene and the camera
% was 370 inch. The physical length of the reference scene was measured by
% placing the tape measure and reading it on the camera captured scene.
% The reference scene was 113.40 inch (h) x 74.80 inch (v).
distanceInch = 370;
refHorizontalInch = 113.40;
refVerticalInch = 74.80;

resolution_camera = [2064 3088];

refVerticalPixel = resolution_camera(1);
refHorizontalPixel = resolution_camera(2);

%% Calculate inch per pixel. We will use this to convert from pixel to inch
% of the target media.
switch options.dir
    case 'vertical'
        pixelToInch = refVerticalInch/refVerticalPixel;
    case 'horizontal'
        pixelToInch = refHorizontalInch/refHorizontalPixel;
end

% Convert the pixel to inch unit.
targetImageInch = pixels_camera * pixelToInch;

%% Calculate the FOV in degrees.
deg = rad2deg(atan(targetImageInch/distanceInch));
fprintf('The FOV of (%d) pixels on the camera is (%.2f) deg. \n',pixels_camera,deg);

%% Calculate the pixel per deg in SACCSFA.
%
% The size of the DMD of the SACCSFA in the camera captured image. The
% resolution of the DMD is [1080 1920].
imageSize_SACCSFA = [1378 2349];
resolution_SACCSFA = [1080 1920];

switch options.dir
    case 'vertical'
        SACCSFAInch = imageSize_SACCSFA(1) * pixelToInch;
    case 'horizontal'
        SACCSFAInch = imageSize_SACCSFA(2) * pixelToInch;
end

deg_SACCSFA = rad2deg(atan(SACCSFAInch/distanceInch));

switch options.dir
    case 'vertical'
        pixelPerDeg_SACCSFA = resolution_SACCSFA(1)/deg_SACCSFA;
    case 'horizontal'
        pixelPerDeg_SACCSFA = resolution_SACCSFA(2)/deg_SACCSFA;
end

pixels_SACCSFA = deg * pixelPerDeg_SACCSFA;
fprintf('On the DMD of the SACCSFA, it is (%.1f) pixels. \n',pixels_SACCSFA);

end
