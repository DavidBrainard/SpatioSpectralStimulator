function [contrasts averagedESF] = GetImgContrast(image, options)
% Calculate the contrast from the given image.
%
% Syntax:
%    contrast = GetImgContrast(image)
%
% Description:
%    It calculates the contrast of input image. The input image should be
%    vertical stripe pattern image.
%
% Inputs:
%    image                      - Image that you want to calculate the
%                                 contrast of.
%
% Outputs:
%    contrast                   - Calculated contrast of the given image.
%
% Optional key/value pairs:
%    verbose                      Boolean. Default true. Controls plotting.
%
% See also:
%    GetImageContrast, GetImageContrastMulti

% History:
%    06/13/23  smo   Wrote it.
%    07/14/23  smo   Updated method to find the valley points on the image.

%% Set variables.
arguments
    image
    options.minPeakDistance (1,1) = 0
    options.verbose (1,1) = true
end

%% Get the image size.
[Ypixel Xpixel] = size(image);

% We will use the average of the 25% / 50% / 75% positions of the cropped image.
imageCrop25 = image(round(0.25*Ypixel),:);
imageCrop50 = image(round(0.50*Ypixel),:);
imageCrop75 = image(round(0.75*Ypixel),:);
imageCropAvg = mean([imageCrop25;imageCrop50;imageCrop75]);

% Find the peaks and valleys. We set the minimum height for both peak and
% valley.
minPeakHeight = max(imageCropAvg) * 0.8;
minValleyHeight = min(-imageCropAvg) * 0.2 + max(-imageCropAvg);
[peaks, locs_peaks] = findpeaks(imageCropAvg,'minpeakdistance',options.minPeakDistance,'minpeakheight',minPeakHeight);
[valleys, locs_valleys] = findpeaks(-imageCropAvg,'minpeakdistance',options.minPeakDistance,'minpeakheight',minValleyHeight);

% Set valley values to the right domain.
valleys = -valleys;

% Calculate contrast here.
nPeaks = length(peaks);
nValleys = length(valleys);

nContrasts = min(nPeaks,nValleys);
for pp = 1:nContrasts-1
    
    % Get contrast per each cycle.
    white = peaks(pp);
    black = valleys(pp);
    contrasts(pp) = (white-black)/(white+black);
end

%% Show the results if you want.
if (options.verbose)
    % Plot it.
    hold on;
    plot(imageCrop25, 'r-', 'LineWidth',1);
    plot(imageCrop50, 'g-', 'LineWidth',1);
    plot(imageCrop75, 'b-', 'LineWidth',1);
    plot(imageCropAvg, 'k-', 'color', [0 0 0 0.2], 'LineWidth',5);
    
    % Plot the peak and valley points.
    plot(locs_peaks,imageCropAvg(locs_peaks),'bo','markersize',4,'markerfacecolor','b','markeredgecolor','k');
    plot(locs_valleys,imageCropAvg(locs_valleys),'ro','markersize',4,'markerfacecolor','r','markeredgecolor','k');
    xlabel('Pixel position (horizontal)','fontsize',15);
    ylabel('dRGB','fontsize',15);
    legend('25%','50%','75%','Avg');
end

% Print out the averaged edge spread function.
averagedESF = imageCropAvg;

end
