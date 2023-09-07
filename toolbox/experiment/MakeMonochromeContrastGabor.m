function [gaborImage, gaborCalUnquantized, gaborCalQuantized, stimulusN, centerN, gaborSizeDeg, gaborSizeMeters] = ...
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
%    gaborCalUnquantized      - The monochrome gabor image in cal format
%                               before quantization.
%    gaborCalQuantized        - The monochrome gabor image in cal format
%                               after quantization. Quantization speeds up
%                               the further image conversion without any
%                               meaningful loss of precision.
%    gaborSizeDeg             - The size of the gabor image in degrees.
%    gaborSizeMeters          - The size of the gabor image in meters.
%    
%
% Optional key/value pairs:
%    nQuantizeBits            - Default to 14. Quantize the contrast image
%                               to a fixed number of levels.
%    sineImagePhaseShiftDeg   - Default to blank. Make a phase shift on the
%                               sine image in degrees. 
%    verbose                  - Deafault true. Print out more status message
%                               if it sets to true.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/20/22  smo             - Wrote it
%   01/25/22  smo             - Added quantization option here.
%   05/09/22  smo             - Added an option to make phase on sine
%                               image.

%% Set parameters.
arguments
    stimulusSizeDeg (1,1)
    sineFreqCyclesPerDeg (1,1)
    gaborSdDeg (1,1)
    screenSizeObj
    options.nQuantizeBits (1,1) = 14
    options.sineImagePhaseShiftDeg = [0]
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

% Make a phase shift on the sine image if you want.
nPhaseShifts = length(options.sineImagePhaseShiftDeg);
if (~isempty(options.sineImagePhaseShiftDeg) & options.sineImagePhaseShiftDeg ~= 0)
    halfStimulusN = stimulusN/2;
    rawMonochromeSineImageSlice = rawMonochromeSineImage(halfStimulusN,:);
    sineAmplitudes = unique(findpeaks(rawMonochromeSineImageSlice));
    
    % Find the peak pixel x-axis coordinates.
    peakLocation = find(ismember(rawMonochromeSineImageSlice,sineAmplitudes));
    oneFrequencyLocation = peakLocation(2) - peakLocation(1);
    
    % Make shift here.
    oneFrequencyDeg = 360;
    amountShift = round(options.sineImagePhaseShiftDeg/oneFrequencyDeg .* oneFrequencyLocation);
    phaseShiftHorizontal = 2;
    for ss = 1:nPhaseShifts
        rawMonochromeSineImagePhaseShift{ss} = circshift(rawMonochromeSineImage, amountShift(ss), phaseShiftHorizontal);
    end
else
    rawMonochromeSineImagePhaseShift{1} = rawMonochromeSineImage;
end

% Make Gaussian window and normalize its max to one.
gaussianWindow = normpdf(MakeRadiusMat(stimulusN,stimulusN,centerN,centerN),0,gaborSdPixels);
gaussianWindow = gaussianWindow/max(gaussianWindow(:));

% Multiply to get Gabor.
for ss = 1:nPhaseShifts
    gaborImage{ss} = rawMonochromeSineImagePhaseShift{ss} .* gaussianWindow;
    
    % Put it into cal format.  Each pixel in cal format is one column. Here
    % there is just one row since it is a monochrome image at this point.
    gaborCalUnquantized{ss} = ImageToCalFormat(gaborImage{ss});
end

if(options.verbose)
    disp('Monochrome gabor image has been created!');
end

%% Quantize the contrast image to a (large) fixed number of levels.
%
% This allows us to speed up the image conversion without any meaningful
% loss of precision. If you don't like it, increase number of quantization
% bits until you are happy again.
nQuantizeLevels = 2^options.nQuantizeBits;
for ss = 1:nPhaseShifts
    gaborCalQuantized{ss} = 2*(PrimariesToIntegerPrimaries((gaborCalUnquantized{ss}+1)/2,nQuantizeLevels)/(nQuantizeLevels-1))-1;
end

if(options.verbose)
    disp('Monochrome gabor cal has been quantized!');
end

% Plot of how quantized version does in obtaining desired contrats.
if (options.verbose)
    for ss = 1:nPhaseShifts
        figure; clf;
        plot(gaborCalUnquantized{ss}(:),gaborCalQuantized{ss}(:),'r+');
        axis('square');
        xlim([0 1]); ylim([0 1]);
        xlabel('Unquantized Gabor contrasts');
        ylabel('Quantized Gabor contrasts');
        title(sprintf('Effect of contrast quantization: Phase shift %d deg\n',options.sineImagePhaseShiftDeg(ss)));
    end
end
end
