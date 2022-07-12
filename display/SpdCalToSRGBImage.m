function [SRGBImage] = SpdCalToSRGBImage(spdImageCal,options)
% Make sRGB image from the spds of image in cal format.
%
% Syntax:
%    [sRGBImage] = SpdCalToSRGBImage(spdImageCal)
%
% Description:
%    This creates sRGB image from the spds of image in cal format. It was
%    written to see the contrast gabor images in standard monitor to check
%    how it looks for SACC project.
%
% Inputs:
%    spdImageCal                - The spds of image in cal format.
%
% Outputs:
%    sRGBImage                  - Created sRGB image in image format.
%
% Optional key/value pairs:
%    S                          - Default to [380 2 201]. The spectral
%                                 range in wavelength for spectrum.
%    verbose                    - Default true. Controls
%                                 plotting and printout.

% History:
%    07/12/22  smo              - Wrote it.

%% Set parameters.
arguments
    spdImageCal
    options.S (1,3) = [380 2 201]
    options.verbose (1,1) = true
end

%% Set image size and load CMFs.
% 
% Set the size of image in pixel. We will make sqaure image.
imageN = sqrt(size(spdImageCal,2));

% Load color matching functions.
load T_xyzJuddVos
T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,options.S);

%% Make SRGB image via XYZ, scaled to display.
predictedXYZCal = T_xyz * spdImageCal;
SRGBPrimaryCal = XYZToSRGBPrimary(predictedXYZCal);
scaleFactor = max(SRGBPrimaryCal(:));
SRGBCal = SRGBGammaCorrect(SRGBPrimaryCal/(2*scaleFactor),0);
SRGBImage = uint8(CalFormatToImage(SRGBCal,imageN,imageN));

%% Show the SRGB image
if (options.verbose)
    figure; imshow(SRGBImage);
    title('SRGB Gabor Image');
end
end
