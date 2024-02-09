function [deg pixels_SACCSFA] = PixelToDeg(pixels_camera,options)
% Convert the unit from pixel to field of view (deg).
%
% Syntax:
%    [deg pixels_SACCSFA] = PixelToDeg(pixels_camera)
%
% Description:
%     It calculates the FOV of the corresponding pixel within the camera
%     captured scene. The calculation is based on the experimental settings
%     that used for the SACC project. Additionally, it calculates the
%     pixel size on the DMD of the SACCSFA system.
%
%     This function was written for quantify the TCA of the SACCSFA system.
%
% Inputs:
%    pixels_camera            - Target length in pixels in the camera
%                               captured scene to convert the unit. It can
%                               be a single number or an array of multiple
%                               numbers.
%
% Outputs:
%    deg                      - Field of view corresponding the input pixel
%                               size in degrees.
%    pixels_SACCSFA           - Corresponding size of the DMD within the
%                               SACCSFA system in pixels.
%
% Optional key/value pairs:
%    dir                      - Default to horizontal. Set the direction of
%                               the screen.
%    verbose                  - Default to true. Controls the printing message.
%
% See also:
%    GetFOVUsingCamera.

% History:
%    01/17/24  smo            - Wrote it.
%    01/18/24  smo            - Updated to work for multiple targets.

%% Set variables.
arguments
    pixels_camera
    options.dir = 'horizontal'
    options.verbose = 'true'
    options.trombone = 'emmetropic'
end

%% Get the reference image info.
%
% The reference image is the image taken far enough to focus with infinity
% camera focus setting. Here, the distance between the scene and the camera
% was 370 inch. The physical length of the reference scene was measured by
% placing the tape measure and reading it on the camera captured scene.
% The reference scene was 113.40 inch (h) x 74.80 inch (v).
distanceInch = 370;
cameraHorizontalInch = 113.40;
cameraVerticalInch = 74.80;

resolution_camera = [2064 3088];

cameraVerticalPixel = resolution_camera(1);
cameraHorizontalPixel = resolution_camera(2);

%% Calculate inch per pixel. We will use this to convert from pixel to inch
% of the target media.
switch options.dir
    case 'vertical'
        pixelToInch = cameraVerticalInch/cameraVerticalPixel;
    case 'horizontal'
        pixelToInch = cameraHorizontalInch/cameraHorizontalPixel;
end

%% Calculate the FOV in degrees.
%
% Convert the pixel to inch unit.
targetImageInch = pixels_camera .* pixelToInch;

% Calculate the FOV in degrees here.
deg = rad2deg(atan(targetImageInch./distanceInch));

%% Calculate the pixel per deg in SACCSFA.
%
% The size of the DMD of the SACCSFA in the camera captured image. The
% resolution of the DMD is [1080 1920].
%
% The DMD size on the camera captured image is different across the
% trombone positions. So, here we set the DMD size on the image differently
% over trombone position. As the trombone position moves from emmetropic
% to compensate near-sighted, the DMD size on the image slightly increases.
switch options.trombone
    case 'emmetropic'
        % This is old measure.
        %         imageSize_SACCSFA = [1378 2349];
        
        % New measure as of 01/24/24.
        imageSize_SACCSFA = [1351 2413];
    case '156'
        imageSize_SACCSFA = [1366 2437];
    case '170'
        imageSize_SACCSFA = [1408 2512];
    case '185'
        imageSize_SACCSFA = [1446 2577];
    otherwise
        error('Trombone position has not been selected properly!');
end
resolution_SACCSFA = [1080 1920];

% Get the FOV of the SACCSFA DMD.
switch options.dir
    case 'vertical'
        SACCSFAInch = imageSize_SACCSFA(1) * pixelToInch;
    case 'horizontal'
        SACCSFAInch = imageSize_SACCSFA(2) * pixelToInch;
end
deg_SACCSFA = rad2deg(atan(SACCSFAInch/distanceInch));

% Calculate the corresponding size in pixel on the SACCSFA DMD.
switch options.dir
    case 'vertical'
        pixelPerDeg_SACCSFA = resolution_SACCSFA(1)/deg_SACCSFA;
    case 'horizontal'
        pixelPerDeg_SACCSFA = resolution_SACCSFA(2)/deg_SACCSFA;
end
pixels_SACCSFA = deg .* pixelPerDeg_SACCSFA;

% Print out the results.
if (options.verbose)
    nTargets = length(deg);
    for tt = 1:nTargets
        % FOV.
        fprintf('The (%s) FOV of (%d) pixels on the camera is (%.2f) deg. \n',...
            options.dir,pixels_camera(tt),deg(tt));
        % Pixels on the DMD.
        fprintf('On the DMD of the SACCSFA, it is (%.1f) pixels. \n\n',pixels_SACCSFA(tt));
    end
end

end
