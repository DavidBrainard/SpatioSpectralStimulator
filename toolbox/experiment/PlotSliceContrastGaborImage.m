function [] = PlotSliceContrastGaborImage(targetGaborImage,desiredGaborImage,options)
% Plot the slice of the center of the gabor contrast image.
%
% Syntax:
%    [] = PlotSliceContrastGaborImage(targetGaborImage,desiredGaborImage)
%
% Description:
%    This plots a slice of the gabor contrast image over the x position
%    (pixels). It takes the target and desired gabor images to compare how
%    well the target image contrast is reproduced as desired image.
%
% Inputs:
%    targetGaborImage              - Target contrast gabor image in image
%                                    format.
%    desiredGaborImage             - Desired contrast gabor image in image
%                                    format. This will be used as reference 
%                                    to compare with the target image.             
%
% Outputs:
%    N/A
%
% Optional key/value pairs:
%    verbose                       - Boolean. Default true. Controls
%                                    plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/26/22  smo                  - Wrote it.

%% Set parameters.
arguments
    targetGaborImage
    desiredGaborImage
    options.plotAxisLimit (1,1) = 4
    options.verbose (1,1) = true
end

%% Plot slice through predicted LMS contrast image.
%
% Note that the y-axis in this plot is individual cone contrast, which is
% not the same as the vector length contrast of the modulation.
%
% Set the size and the center point of the gabor image.
if (any(size(targetGaborImage) ~= size(desiredGaborImage))) 
    error('Target and desired gabor image size do not match!');
end
stimulusN = size(targetGaborImage,1);
centerN = stimulusN * 0.5;

% Plot it.
figure; clf; hold on
plot(1:stimulusN,100 * targetGaborImage(centerN,:,1),'ro','MarkerFaceColor','r','MarkerSize',2);
plot(1:stimulusN,100 * desiredGaborImage(centerN,:,1),'r','LineWidth',0.5);

plot(1:stimulusN,100 * targetGaborImage(centerN,:,2),'go','MarkerFaceColor','g','MarkerSize',2);
plot(1:stimulusN,100 * desiredGaborImage(centerN,:,2),'g','LineWidth',0.5);

plot(1:stimulusN,100 * targetGaborImage(centerN,:,3),'bo','MarkerFaceColor','b','MarkerSize',2);
plot(1:stimulusN,100 * desiredGaborImage(centerN,:,3),'b','LineWidth',0.5);
ylabel('LMS Cone Contrast (%)');

if (size(targetGaborImage,3) == 4)
    plot(1:stimulusN,100 * targetGaborImage(centerN,:,4),'co','MarkerFaceColor','b','MarkerSize',2);
    plot(1:stimulusN,100 * desiredGaborImage(centerN,:,4),'c','LineWidth',0.5);
    ylabel('LMS Cone and Mel Contrast (%)');
end

xlabel('x position (pixels)')
ylim([-options.plotAxisLimit options.plotAxisLimit]);

end