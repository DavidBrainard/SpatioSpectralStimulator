function [imageConverted rngVal] = AddNoiseToImage(image, options)
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
%    imageConverted           - Output image that added noise to the given
%                               image.
%    rngVal                   - Seed of the random noise using rng fucntion
%                               so that we can recreate the noise pattern
%                               when we want.
%
% Optional key/value pairs:
%    noiseLevel               - Default to 1. The maximum noise intensity
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
%    08/12/22  smo            - Started on it.
%    08/17/22  smo            - Updated on sampling. Now we can make noise
%                               image for any noise intensity.
%    08/18/22  dhb, smo       - Now we add the noise without the specific
%                               pattern to all pixels.
%    09/08/22  smo            - We save the seed number for random noise.

%% Set parameters.
arguments
    image
    options.noiseLevel (1,1) = 1
    options.nInputLevels (1,1) = 256
    options.verbose (1,1) = false
end

%% Make noise image, add to input, and handle out of range.
%
% Save the seed for random noise so that we can recreate it if we want.
rngVal = rng;

% Make a noise image.
noiseImage = randi([-options.noiseLevel options.noiseLevel], size(image));

% Add noise to image.
imageConverted = image + noiseImage./options.nInputLevels;

% Handle out of range.
imageConverted(imageConverted < 0) = 0;
imageConverted(imageConverted > 1) = 1;

%% OLD LINES.
% % Set the noise intensity parameters.
% intervalNoiseIntensity = 1;
% intensityNoise = [-options.noiseLevel: intervalNoiseIntensity : options.noiseLevel]./options.nInputLevels;
% nIntensityNoise = length(intensityNoise);
% numelRangeNoise = round(imageN/nIntensityNoise);
% 
% % Separate the ranges randomly.
% rangeNoise = 1:imageN;
% rangeNoiseTemp = [];
% for nn = 1:nIntensityNoise
%     noiseSamplingRange = setdiff(rangeNoise,rangeNoiseTemp);
%     
%     % If there are not enough elements left for sampling, just fill out
%     % with the first pixels.
%     if (numelRangeNoise > length(noiseSamplingRange))
%         numElementsToAdd = numelRangeNoise - length(noiseSamplingRange) + 1;
%         noiseSamplingRange = [noiseSamplingRange ones(1, numElementsToAdd)];
%     end
%     
%     rangeNoiseSet(:,nn) = randsample(noiseSamplingRange, numelRangeNoise);
%     rangeNoiseTemp = rangeNoiseSet(:);
% end
% 
% % Color the noise image.
% for nn = 1:nIntensityNoise
%     noiseImage(rangeNoiseSet(:,nn), :, whichPrimary) = intensityNoise(nn);
% end
% 
% %% Add noise to the given image.
% imageConverted = image + noiseImage;

%% Show the image added noise.
if (options.verbose)
    figure;
    imshow(imageConverted);
    title(sprintf('Noise Intensity = %d', options.noiseLevel), 'FontSize', 15);
end

end
