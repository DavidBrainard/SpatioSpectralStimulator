function contrast = GetImgContrast(image, options)
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
    options.verbose (1,1) = true
end

%% Get the image size.
[Ypixel Xpixel] = size(image);

% We will use the average of the 25% / 50% / 75% positions of the cropped image.
imageCrop25 = image(round(0.25*Ypixel),:);
imageCrop50 = image(round(0.50*Ypixel),:);
imageCrop75 = image(round(0.75*Ypixel),:);
imageCropAvg = mean([imageCrop25;imageCrop50;imageCrop75]);

% Calculate contrast here.
white = max(imageCropAvg);
black = min(imageCropAvg);
contrast = (white-black)/(white+black);

%% Show the results if you want.
if (options.verbose)
    % Show image.
%     figure; imshow(image);
    
    % Plot it.
%     figure; hold on;
      hold on;
    plot(imageCrop25, 'r-', 'LineWidth',1);
    plot(imageCrop50, 'g-', 'LineWidth',1);
    plot(imageCrop75, 'b-', 'LineWidth',1);
    plot(imageCropAvg, 'k-', 'color', [0 0 0 0.2], 'LineWidth',5);
%     plot(imageCropAvg, 'k-', 'LineWidth',1);
    xlabel('Pixel position (horizontal)','fontsize',15);
    ylabel('dRGB','fontsize',15);
    legend('25%','50%','75%','Avg');
%     legend(sprintf('Contrast = %.2f',contrast));
end
end
