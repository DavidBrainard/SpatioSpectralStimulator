% GetImageContrastMulti
%
% This is to calculate the contrast of the image taken by camera. This is
% to calculate multiple images at once.

% History:
%    09/16/22  smo   Wrote it.

%% Initialize.
clear; close all; 

%% Load image.
% 
% Direct to the folder.
fileDir = '~/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Camera/091622_Check_Chromatic_Abberration';
cd(fileDir);

% Load image here.
for ii = 1:16
    fileName = sprintf('ch%d_FocusFixed',ii);
    imageFileName = GetMostRecentFileName(fileDir,fileName);
    image{ii} = imread(imageFileName);
end 

%% Find the rectangle on the image.
FINDRECTANGLE = false;

if (FINDRECTANGLE)
    % Get the size of the image.
    [Ypixel Xpixel] = size(image); 
    
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
    
    % Show the image and draw the detected rectangles on it
    figure; imshow(image); hold on;
    for i = 1:numel(stats)
        rectangle('Position',stats(i).BoundingBox,'Linewidth',3,'EdgeColor','y','LineStyle','-');
        area(i) = stats(i).Area;
        centroid(i,:) = stats(i).Centroid;
    end
    
    % Spot the target on the image.
    targetStat = 6;
    rectangle('Position',stats(targetStat).BoundingBox,'Linewidth',3,'EdgeColor','r','LineStyle','--');
end

%% Clip the part of contrast image.
fromX = 1522;
fromY = 1096;
cropXPixel = 70;
cropYPixel = 31;
toX = fromX+70;
toY = fromY+31;

for ii = 1:16
    imageTemp = image{ii};
    imageCropTemp = imageTemp(fromY:toY,fromX:toX);
    [YpixelCrop XpixelCrop] = size(imageCropTemp);
    
    % We will use the average of the 25% / 50% / 75% positions of the cropped image.
    imageCrop25(:,ii) = imageCropTemp(round(0.25*YpixelCrop),:);
    imageCrop50(:,ii) = imageCropTemp(round(0.50*YpixelCrop),:);
    imageCrop75(:,ii) = imageCropTemp(round(0.75*YpixelCrop),:);
    imageCropAvg = mean([imageCrop25;imageCrop50;imageCrop75]);
    
    % Calculate the contrast here.
    whiteCropImage = max(im2double(imageCrop50(:,ii)));
    blackCropImage = min(im2double(imageCrop50(:,ii)));
    contrast(ii) = (whiteCropImage-blackCropImage) / (whiteCropImage+blackCropImage);
end

%% Plot it.
figure; hold on;
plot(imageCrop50, 'LineWidth',1);
xlabel('Pixel position (horizontal)','fontsize',15);
ylabel('dRGB','fontsize',15);
leg = string([1:16]);
legend(leg,'fontsize',12);
