function [imageConverted] = AddNoiseToImage(image, options)
% Add noise to given image.
%
% Syntax:
%    [noiseImage] = AddNoiseToImage(image)
%
% Description:
%    This adds noise to image. As we observe artifacts on the contrast
%    images in the SACC project, we try to add random noise to see if that
%    makes it to perceive less artifacts. 
%
%    Strategy here is to add random noise within the set number of maximum
%    intensity (default to 1), and allocate the noise with different
%    intensities on the image at randomly sampled range. So, the noise
%    itself would look like un-even stripes. 
%
% Inputs:
%    image                    - Image that you want to add noise to.
%
% Outputs:
%    noiseImage               - Output image that added noise to the given
%                               image.
%
% Optional key/value pairs:
%    noiseIntensityMax        - Default to 1. The maximum noise intensity
%                               in dRGB unit. The random noise will be
%                               generated in the range within
%                               [-this value : 1 : this value]. So, if it
%                               set to 1, the range would be [-1, 0, 1].
%    nInputLevels             - Default to 256. The number of input levels
%                               on the system. Default value follows the
%                               SACC projector.
%    verbose                  - Default to false. Controls plotting.
%
% See also:
%    N/A

% History:
%   08/12/22  smo             - Started on it.

%% Set parameters.
arguments
    image
    options.noiseIntensityMax (1,1) = 2
    options.nInputLevels (1,1) = 256
    options.verbose (1,1) = false
end

%% Get the size of the image. 
imageN = size(image,1);
nPrimaries = size(image,3);

%% Make noise image.
noiseImage = zeros(imageN, imageN, nPrimaries);
whichPrimary = 1:nPrimaries;

% Set the noise intensity parameters.
intervalNoiseIntensity = 1;
intensityNoise = [-options.noiseIntensityMax: intervalNoiseIntensity : options.noiseIntensityMax]./options.nInputLevels;
nIntensityNoise = length(intensityNoise);
numelRangeNoise = round(imageN/nIntensityNoise);

% Separate the ranges randomly.
rangeNoise = 1:imageN;
rangeNoiseTemp = [];
for nn = 1:nIntensityNoise    
    rangeNoiseSet(:,nn) = randsample(setdiff(rangeNoise,rangeNoiseTemp), numelRangeNoise);
    rangeNoiseTemp = rangeNoiseSet(:);
end

% Color the noise image.
for nn = 1:nIntensityNoise
    noiseImage(rangeNoiseSet(:,nn), :, whichPrimary) = intensityNoise(nn);
end

%% Add noise to the given image.
imageConverted = image + noiseImage;

% Show the image added noise.
if (options.verbose)
    figure;
    imshow(imageConverted);
    title(sprintf('Noise Intensity = %d', options.noiseIntensityMax), 'FontSize', 15);
end

end
