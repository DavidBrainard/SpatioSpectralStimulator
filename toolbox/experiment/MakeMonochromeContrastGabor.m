function [gaborImage, gaborCal, stimulusN, centerN, gaborSizeDeg, gaborSizeMeters] = ...
    MakeMonochromeContrastGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenSizeObj,options)
% Make a monochrome gabor patch image.
%
% Syntax:
%    [gaborImage, gaborCal, stimulusN, centerN, gaborSizeDeg, gaborSizeMeters] = ...
%       MakeMonochromeGabor(stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenPixelsPerDeg,screenDpm)
%
% Description:
%    This makes a monochrome gabor patch image. The shape of gabor image is
%    decided by input parameters and this is pre-step of making the gabor
%    patch modulated in a color direction.
%
% Inputs:
%    stimulusSizeDeg          - Size of the gabor image in degrees.
%    sineFreqCyclesPerDeg     - The number of sine cycles per each degree
%                               in the gabor image.
%    gaborSdDeg               - Standard deviation of the gabor image in degrees.
%    screenPixelsPerDeg       - The number of screen pixels per degrees.
%    screenDpm                - The number of screen dots per meters (DPM).
%
% Outputs:
%    gaborImage               - Created monochrome gabor image in image
%                               format.
%    gaborCal                 - The monochrome gabor image in cal format.
%    gaborSizeDeg             - The size of the gabor image in degrees.
%    gaborSizeMeters          - The size of the gabor image in meters.
%
% Optional key/value pairs:
%    verbose                  - Deafault true. Print out more status message
%                               if it sets to true.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/20/22  smo             - Wrote it

%% Set parameters.
arguments
    stimulusSizeDeg (1,1)
    sineFreqCyclesPerDeg (1,1)
    gaborSdDeg (1,1)
    screenSizeObj
    options.verbose (1,1) = true
end

%% Say hello.
if(options.verbose)
    fprintf('Making Gabor contrast image\n');
end

%% Set stimulus size in pixels.
%
% Stimulus goes into a square image.  Want number of pixels to be even. If
% we adjust pixels, also adjust image size in degrees.
stimulusN = round(vecnorm([stimulusSizeDeg stimulusSizeDeg])*screenSizeObj.screenPixelsPerDeg/sqrt(2));
if (rem(stimulusN,2) ~= 0)
    stimulusN = stimulusN+1;
end
centerN = stimulusN/2;

% Compute expected stimulus size in degrees based on actual pixels in the
% square image.
gaborSizeDeg = vecnorm([stimulusN stimulusN])/(screenSizeObj.screenPixelsPerDeg*sqrt(2));
gaborSizeMeters = stimulusN/screenSizeObj.screenDpm;

%% Convert image parameters in degrees to those in pixels.
sineFreqCyclesPerImage = sineFreqCyclesPerDeg * stimulusSizeDeg;
gaborSdPixels = gaborSdDeg * screenSizeObj.screenPixelsPerDeg;

%% Make sine image here.
%
% The function MakeSineImage actually makes plaids. By making the
% horozontal frequency (first argument) 0, we get a vertical sinusoid.
rawMonochromeSineImage = MakeSineImage(0,sineFreqCyclesPerImage,stimulusN);

% Make Gaussian window and normalize its max to one.
gaussianWindow = normpdf(MakeRadiusMat(stimulusN,stimulusN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));

% Multiply to get Gabor.
gaborImage = rawMonochromeSineImage.*gaussianWindow;

% Put it into cal format.  Each pixel in cal format is one column. Here
% there is just one row since it is a monochrome image at this point.
gaborCal = ImageToCalFormat(gaborImage);

if(options.verbose)
    disp('Monochrome gabor image has been created!');
end
end
