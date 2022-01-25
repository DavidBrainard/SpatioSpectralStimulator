function [ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj,options)
% Set up ISETBio display object from the screen cal object.
%
% Syntax:
%    [ISETBioDisplayObject,screenSizeObject,screenCalObjFromISETBio] = SetupISETBioDisplayObject(colorDirectionParams,screenCalObj)
%
% Description:
%    It creates ISETBio display object from screen cal object as we want to
%    utilize ISETBio in the SACC project. It also prints out the screen
%    calibration object as an output, which is taken back from ISETBio
%    object to make sure that we can move the domain from one to another
%    with no problems.
%
% Inputs:
%    colorDirectionParams      - Structure with the parameters to
%                                calculate a contrast gabor image.
%    screenCalObj              - Screen calibration object.
%
% Outputs:
%    ISETBioDisplayObject      - Structure with the parameters to make the
%                                ISETBio scene from image.
%    screenSizeObject          - Structure that contains the target display
%                                size in different units such as meters,
%                                inches, degs, and pixels.
%    screenCalObjFromISETBio   - Screen calibration object that is
%                                calculated from the ISETBio display
%                                object. This should be matched with the
%                                initial screen cal object.
%
% Optional key/value pairs:
%    verbose                   - Boolean. Default true. Controls
%                                plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,ga,smo     - Wrote it.
%   01/24/22  smo            - Made it work.

%% Set parameters.
arguments
    colorDirectionParams
    screenCalObj
    options.verbose (1,1) = true
end

%% Set the display size in different units.
%
% Set the display parameters.
screenDiagSizeDeg = 15.5;
screenDistanceMeters = 10;
screenSizePixels = [1920 1080];
screenDiagSizePixels = vecnorm(screenSizePixels);

% Convert the screen size from degrees to meters/inches.
%
% Note that we need to do this along the diagonal because degrees aren't
% linear in meters, so we want to work first in the physical units of the
% display, not in degrees.
[screenDiagSizeMeters,screenSizeMeters,screenSizeInches] = ...
    DegToMeters(screenDiagSizeDeg,screenDistanceMeters,screenSizePixels,'DegToInches',true);

% Get horizontal and vertical size of screen in degrees. We take pixels per
% degree along the diagonal as the best compromise, need to use that when
% we compute image sizes below.
screenSizeDeg = 2 * atand(screenSizeMeters/(2*screenDistanceMeters));
screenPixelsPerDeg = screenDiagSizePixels / screenDiagSizeDeg;

% Get dpi and make sure everything is consistent.
screenDpi = vecnorm(screenSizePixels)/vecnorm(screenSizeInches);
screenDpiChk = mean(screenSizePixels ./ screenSizeInches);
if (abs((screenDpi - vecnorm(screenDpiChk))/screenDpi) > 1e-6)
    error('Screen is not rigid');
end
inchesPerMeter = 39.3701;
screenDpm = screenDpi * inchesPerMeter;

%% Create ISETBio display object here.
%
% Note that the function takes T_cones, S, screenGammaMethod for checking
% if the reverse calculation (from ISETBio to calibration object) can get
% the same values.
screenGammaMethod = 2;
[ISETBioDisplayObject,screenCalObjFromISETBio] = MakeISETBioDisplayObj(screenCalObj,screenDistanceMeters,...
    screenSizeMeters,screenSizePixels,colorDirectionParams.T_cones,colorDirectionParams.S,screenGammaMethod,'verbose',options.verbose);

%% Save the screen size params in a struct.
screenSizeObject.screenDiagSizeDeg = screenDiagSizeDeg;
screenSizeObject.screenDistanceMeters = screenDistanceMeters;
screenSizeObject.screenSizeMeters = screenSizeMeters;
screenSizeObject.screenDiagSizeMeters = screenDiagSizeMeters;
screenSizeObject.screenSizeInches = screenSizeInches;
screenSizeObject.screenSizePixels = screenSizePixels;
screenSizeObject.screenDiagSizePixels = screenDiagSizePixels;
screenSizeObject.screenSizeDeg = screenSizeDeg;
screenSizeObject.screenPixelsPerDeg = screenPixelsPerDeg;
screenSizeObject.screenDpi = screenDpi;
screenSizeObject.screenDpm = screenDpm;

end
