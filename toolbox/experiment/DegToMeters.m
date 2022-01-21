function [screenDiagSizeMeters, screenSizeMeters, screenSizeInches] = DegToMeters(screenDiagSizeDeg,screenDistanceMeters,screenSizePixels,options)
% Convert the screen size from degree to meters.
%
% Syntax:
%    [screenDiagSizeMeters, screenSizeMeters, screenSizeInches] = DegToMeters(screenDiagSizeDeg,screenDistanceMeters,screenSizePixels)
%
% Description:
%    This is for converting the screen size from visaul angle (deg) to
%    physical unit (meters, inches). This is especially for the SACC
%    project using the ISETBio, which takes the screen size information in
%    meters. 
%
% Inputs:
%    screenDiagSizeDeg        - Target screen diagonal size in degrees.
%    screenDistanceMeters     - Distance between the screen and the
%                               observer. For SACC project, it may set as
%                               high enough number to represent the
%                               infinity.
%    screenSizePixels         - Screen resolution in pixels. If the
%                               resolution is 1920 x 1080, it should be
%                               typed as [1920 1080].
%
% Outputs:
%    screenDiagSizeMeters     - Target screen diagonal size in meters.
%    screenSizeMeters         - Screen size in meters as [horizontal vertical].
%    screenSizeInches         - Screen size in inches as [horizontal vertical]. 
%
% Optional key/value pairs:
%    DegToInches              - Default true. If it is set as true, you
%                               will also get the results of the screen
%                               size in inches.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/20/22  smo             - Wrote it

%% Set parameters.
arguments
    screenDiagSizeDeg (1,1)
    screenDistanceMeters (1,1)
    screenSizePixels (1,2)
    options.DegToInches (1,1) = true
end

%% Get screen pixel values for horizontal and vertical sizes.
screenHorizSizePixels = screenSizePixels(1);
screenVertSizePixels  = screenSizePixels(2);

%% Convert the size here.
screenDiagSizeMeters = 2 * screenDistanceMeters * tand(screenDiagSizeDeg/2);
screenHorizSizeMeters = sqrt(screenDiagSizeMeters^2 / (1+(screenVertSizePixels/screenHorizSizePixels)^2));
screenVertSizeMeters = screenHorizSizeMeters * screenVertSizePixels/screenHorizSizePixels;
screenSizeMeters = [screenHorizSizeMeters screenVertSizeMeters];

% Check the calculation.
if (abs((screenDiagSizeMeters - vecnorm(screenSizeMeters))/screenDiagSizeMeters) > 1e-6)
    error('You did not understand what Pythagoras said!');
end

%% Calculate the meters to inches if you want.
inchesPerMeter = 39.3701;
if (options.DegToInches)
    screenSizeInches = screenSizeMeters * inchesPerMeter;
end

end