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
    options.xAxisValues = [];
    options.xAxisLabel = [];
    options.markerSize = 2;
    options.lineWidth = 0.5;
    options.plotAxisPos (1,1) = false;
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

% Set up x-axis values
if isempty(options.xAxisValues)
    xAxisValues = 1:stimulusN;
else
    xAxisValues = options.xAxisValues;
end
if (isempty(options.xAxisLabel))
    xAxisLabel = 'Position (pixels)';
else
    xAxisLabel = options.xAxisLabel;
end

% Plot it.
figure; clf; hold on
plot(xAxisValues,100 * desiredGaborImage(centerN,:,1),'r','LineWidth',options.lineWidth);
plot(xAxisValues,100 * desiredGaborImage(centerN,:,2),'g','LineWidth',options.lineWidth);
plot(xAxisValues,100 * desiredGaborImage(centerN,:,3),'b','LineWidth',options.lineWidth);
if (size(targetGaborImage,3) == 4)
    plot(xAxisValues,100 * desiredGaborImage(centerN,:,4),'c','LineWidth',options.lineWidth);
    ylabel('LMS Cone and Mel Contrast (%)');
end

if (options.markerSize ~= 0)
    plot(xAxisValues,100 * targetGaborImage(centerN,:,1),'ro','MarkerFaceColor','r','MarkerSize',options.markerSize);
    plot(xAxisValues,100 * targetGaborImage(centerN,:,2),'go','MarkerFaceColor','g','MarkerSize',options.markerSize);
    plot(xAxisValues,100 * targetGaborImage(centerN,:,3),'bo','MarkerFaceColor','b','MarkerSize',options.markerSize);
    if (size(targetGaborImage,3) == 4)
        plot(xAxisValues,100 * targetGaborImage(centerN,:,4),'co','MarkerFaceColor','b','MarkerSize',options.markerSize);
    end
end

if (size(targetGaborImage,3) == 4)
    ylabel('LMS Cone and Mel Contrast (%)');
else
    ylabel('LMS Cone Contrast (%)');
end



xlabel(xAxisLabel)
if (options.plotAxisPos)
    ylim([0 options.plotAxisLimit]);
else
    ylim([-options.plotAxisLimit options.plotAxisLimit]);
end

end