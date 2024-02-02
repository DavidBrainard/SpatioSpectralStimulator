function [contrast] = GetIPContrast(IP, options)
% Calculate the contrast from the intensity profile across horizontal pixel positions.
%
% Syntax:
%    contrast = GetImgContrast(image)
%
% Description:
%    It calculates the contrast of input image. The input image should be
%    vertical stripe pattern image.
%
% Inputs:
%    IP                         - Target image's intensity profile across
%                                 horizontal pixels that you want to
%                                 calculate the contrast of.
%
% Outputs:
%    contrast                   - Calculated contrast of the given
%                                 intensity profile of the image.
%
% Optional key/value pairs:
%    verbose                      Boolean. Default true. Controls plotting.
%
% See also:
%    GetImgContrast, GetImageContrastMulti

% History:
%    01/03/24    smo            - Wrote it. This is basically the same
%                                 function as 'GetImgContrast', but it
%                                 takes intensity profile rather than the
%                                 raw image.

%% Set variables.
arguments
    IP
    options.minPeakDistance (1,1) = 0
    options.verbose (1,1) = true
end

%% Check the class of the IP and switch it to double if it's not. Further calculations should be done in double.
IP_class = class(IP);

if ~strcmp(IP_class,'double')
    IP = double(IP);
end

%% Find the peaks and valleys. We set the minimum height for both peak and
% valley.
minPeakHeight = max(IP) * 0.8;
minValleyHeight = min(-IP) * 0.2 + max(-IP);
[peaks, locs_peaks] = findpeaks(IP,'minpeakdistance',options.minPeakDistance,'minpeakheight',minPeakHeight);
[valleys, locs_valleys] = findpeaks(-IP,'minpeakdistance',options.minPeakDistance,'minpeakheight',minValleyHeight);

% Set valley values to the right domain.
valleys = -valleys;

% Get the number of peaks and vallyes in the given range..
nPeaks = length(peaks);
nValleys = length(valleys);

% Based on the number of peak and valley, we set the number of contrasts to
% calculate. We follow whichever a smaller number between peak and valley.
nContrasts = min(nPeaks,nValleys);

% Calculate contrast here.
for pp = 1:nContrasts
    
    % Get contrast per each cycle.
    white = peaks(pp);
    black = valleys(pp);
    contrastsRaw(pp) = (white-black)/(white+black);
end

% Calculate the mean contrasts.
contrast = mean(contrastsRaw);

%% Show the results if you want.
if (options.verbose)
    % Plot it.
    plot(IP, 'k-', 'color', [0 0 0 0.2], 'LineWidth',5);
    
    % Plot the peak and valley points.
    plot(locs_peaks,IP(locs_peaks),'bo','markersize',4,'markerfacecolor','b','markeredgecolor','k');
    plot(locs_valleys,IP(locs_valleys),'ro','markersize',4,'markerfacecolor','r','markeredgecolor','k');
    xlabel('Pixel position (horizontal)','fontsize',15);
    ylabel('dRGB','fontsize',15);
    legend('Intensity profile');
end

end
