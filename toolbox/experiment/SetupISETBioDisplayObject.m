function [] = SetupISETBioDisplayObject()
% d
%
% Syntax:
%    d
%
% Description:
%    d
%
% Inputs:
%    d                       -
%
% Outputs:
%    d                       -
%
% Optional key/value pairs:
%    d                       - d
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,ga,smo     - Wrote it

%% Set parameters.
arguments
end

%%
% Note that we will want to be able to pass in custom primaries to this
% routine, or at least do it with different cal files, because we'll do it
% sometimes with nominal primaries and sometimes with measured primaries.

%
% On the DMD size, Derek writes (email 2022-01-11):
%     Looking at the system design, I find that the conversion is 7.74350
%     degrees : 12 mm = 0.64529 degrees per mm. This indicates that the
%     diagonal of the DMD occupies a FOV of 2*7.74350 = 15.487 degrees ~
%     15.5 degrees.
%
%     We can use this to solve for the h and v size in degrees.
%     Let the horizontal size by x and the diagonal be d.  We know
%     that d^2 = x^2 * (1+(screenVertPixels/screenHorizPixels)^2). So
%     x = sqrt(d^2/(1+(screenVertPixels/screenHorizPixels)^2)).
%
% The DMD dimensions and distance are dummied up so that the visual
% angle matches that of our optical system, but so that the optical
% distance is large enough to mimic optical infinity.
%
% Display parameters.
screenDiagSizeDeg = 15.5;
screenDistanceVirtualMeters = 10;
screenSizePixels = [1920 1080];
screenDiagSizePixels = vecnorm(screenSizePixels);

% Convert the screen size from degrees to meters/inches.
%
% Note that we need to do this along the diagonal because degrees aren't
% linear in meters, so we want to work first in the physical units of the
% display, not in degrees.
[screenDiagSizeMeters,screenSizeMeters,screenSizeInches] = ...
    DegToMeters(screenDiagSizeDeg,screenDistanceVirtualMeters,screenSizePixels,'DegToInches',true);

% Get horizontal and vertical size of screen in degrees. We take pixels per
% degree along the diagonal as the best compromise, need to use that when
% we compute image sizes below.
screenSizeDeg = 2 * atand(screenSizeMeters/(2*screenDistanceVirtualMeters));
screenPixelsPerDeg = screenDiagSizePixels / screenDiagSizeDeg;

% Get dpi and make sure everything is consistent.
screenDpi = vecnorm(screenSizePixels)/vecnorm(screenSizeInches);
screenDpiChk = mean(screenSizePixels ./ screenSizeInches);
if (abs((screenDpi - vecnorm(screenDpiChk))/screenDpi) > 1e-6)
    error('Screen is not rigid');
end
inchesPerMeter = 39.3701;
screenDpm = screenDpi*inchesPerMeter;

% Create ISETBio display object here.
%
% Note that the function takes T_cones, S, screenGammaMethod for checking
% if the reverse calculation (from ISETBio to calibration object) can get
% the same values.
[ISETBioDisplayObject,screenCalObjFromISETBio] = MakeISETBioDisplayObj(screenCalObj,screenDistanceVirtualMeters,...
    screenSizeMeters,screenSizePixels,colorDirectionParams.T_cones,colorDirectionParams.S,screenGammaMethod,'verbose',VERBOSE);
end
