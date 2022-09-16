% GetImageContrast
%
% This is to calculate the contrast of the image taken by camera.

% History:
%    09/16/22  smo   Wrote it.

%% Initialize.
clear; close all; 

%% Load image.
fileDir = '~/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Camera/091622_Check_Chromatic_Abberration';
cd(fileDir);

numChannel = 16;
fileName = sprintf('ch%d_FocusNot',numChannel);
imageFileName = GetMostRecentFileName(fileDir,fileName);
image = imread(imageFileName);
[Ypixel Xpixel] = size(image); 

%% Find the rectangle on the image.
FINDRECTANGLE = false;

if (FINDRECTANGLE)
    % find both black and white regions
    stats = [regionprops(image); regionprops(not(image))]; % Find the objects on the image
    CentroidDouble = cell2mat({stats.Centroid}.');
    
    % Centroid limit
    lowerlimitXpixel = 0.46;
    upperlimitXpixel = 0.54;
    lowerlimitYpixel = 0.46;
    upperlimitYpixel = 0.54;
    indexcentroid = find(CentroidDouble(:,1) > Xpixel*lowerlimitXpixel & CentroidDouble(:,1) < Xpixel*upperlimitXpixel & CentroidDouble(:,2) > Ypixel*lowerlimitYpixel & CentroidDouble(:,2) < Ypixel*upperlimitYpixel);
    stats = stats(indexcentroid);
    
    % Area limit
    AreaDouble = cell2mat({stats.Area}.');
    lowerlimitarea = 40;
    upperlimitarea = 60;
    indexarea = find(AreaDouble > lowerlimitarea & AreaDouble < upperlimitarea);
    stats = stats(indexarea);
    
    % show the image and draw the detected rectangles on it
    figure; imshow(image); hold on;
    for i = 1:numel(stats)
        rectangle('Position',stats(i).BoundingBox,'Linewidth',3,'EdgeColor','y','LineStyle','-');
        area(i) = stats(i).Area;
        centroid(i,:) = stats(i).Centroid;
    end
    
    targetStat = 6;
    rectangle('Position',stats(targetStat).BoundingBox,'Linewidth',3,'EdgeColor','r','LineStyle','--');
end

%% Contrast measurement routine.
imageSize = size(image);
screenXpixel = imageSize(1);
screenYpixel = imageSize(2);

% Clip the part of the screen image to calculate the contrast.
fromX = 1522;
fromY = 1096;
cropXPixel = 70;
cropYPixel = 31;
toX = fromX+70;
toY = fromY+31;
imageCrop = image(fromY:toY,fromX:toX);
[YpixelCrop XpixelCrop] = size(imageCrop);

% We will use the average of the 25% / 50% / 75% positions of the cropped image. 
imageCrop25 = imageCrop(round(0.25*YpixelCrop),:);
imageCrop50 = imageCrop(round(0.50*YpixelCrop),:);
imageCrop75 = imageCrop(round(0.75*YpixelCrop),:);
imageCropAvg = mean([imageCrop25;imageCrop50;imageCrop75]);

% Plot it.
figure; hold on;
plot(imageCrop25, 'r-', 'LineWidth',1);
plot(imageCrop50, 'g-', 'LineWidth',1);
plot(imageCrop75, 'b-', 'LineWidth',1);
plot(imageCropAvg, 'k--', 'LineWidth',2.5);
xlabel('Pixel position (horizontal)');
ylabel('dRGB');
legend('25%','50%','75%','Avg');

% Calculate contrast here.
whiteCropImage = max(imageCropAvg);
blackCropImage = min(imageCropAvg);
contrastCropImage = (whiteCropImage-blackCropImage) / (whiteCropImage+blackCropImage);

contrastCropImage
