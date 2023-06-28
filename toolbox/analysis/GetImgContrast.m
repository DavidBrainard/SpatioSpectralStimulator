function contrasts = GetImgContrast(image, options)
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

% Find the peaks.
minPeakHeight = max(imageCropAvg) * 0.8;
findpeaks(imageCropAvg,'minpeakdistance',options.minPeakDistance,'minpeakheight',minPeakHeight);

% Calculate contrast here.
%
% We will calculate the contrast per each cycle.
[peaks, locs_peaks] = findpeaks(imageCropAvg,'minpeakdistance',options.minPeakDistance,'minpeakheight',minPeakHeight);

nPeaks = length(peaks);
for pp = 1:nPeaks-1
    % Find the valley locations.
    locs_valleys(pp) = round(locs_peaks(pp) + (locs_peaks(pp+1)-locs_peaks(pp))/2);

    % Calculate contrast here.
    white = peaks(pp);
    black = imageCropAvg(locs_valleys(pp));
    contrasts(pp) = (white-black)/(white+black);
end

% white = max(imageCropAvg);
% black = min(imageCropAvg);
% contrast = (white-black)/(white+black);

%% Show the results if you want.
if (options.verbose)
    % Show image.
%     imshow(image);

    % Plot it.
    hold on;
    plot(imageCrop25, 'r-', 'LineWidth',1);
    plot(imageCrop50, 'g-', 'LineWidth',1);
    plot(imageCrop75, 'b-', 'LineWidth',1);
    plot(imageCropAvg, 'k-', 'color', [0 0 0 0.2], 'LineWidth',5);

    % Plot the valley points.
    plot(locs_valleys,imageCropAvg(locs_valleys),'ro','markersize',4,'markerfacecolor','r')
    xlabel('Pixel position (horizontal)','fontsize',15);
    ylabel('dRGB','fontsize',15);
    legend('25%','50%','75%','Avg');
end
end
