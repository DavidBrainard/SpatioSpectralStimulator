%% Initialize
clear; close all;

%% Load an image and find targeted circle
image = imread('Test3.tiff','tif');
[Ypixel Xpixel] = size(image);

figure(1); imshow(image); title('Original Image');

% Find size of circles we want
d = imdistline; % measure the radius of the target circle image

image_filtered = imfilter(image,ones(5)/25);

% Test 1 center image circle
% [centers, radii] = imfindcircles(image_filtered,[90 130]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
% Test 1 big circle
% [centers, radii] = imfindcircles(image_filtered,[1400 1500]/2,'ObjectPolarity','dark','Sensitivity',0.98,'EdgeThreshold',0.03);
% Test 3 image
[centers, radii] = imfindcircles(image_filtered,[120 160]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);

% LCPA1 image
% [centers, radii] = imfindcircles(image_filtered,[10 20]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
% [10 18] for the position 1
% [14 22] for the position 2

h = viscircles(centers,radii);

% figure; imshow(image);
% [centers2, radii2] = imfindcircles(image_filtered,[40 60]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
% h2 = viscircles(centers2,radii2);

% Find targeted circles 
figure(2); imshow(image_filtered); title('Image with targeted points'); hold on;

% % Centered circle on the test pattern
lowerlimit_X = 0.1;
upperlimit_X = 0.8;
lowerlimit_Y = 0.1;
upperlimit_Y = 0.9;

% % LCPA1 centered point 
% lowerlimit_X = 0.45;
% upperlimit_X = 0.52;
% lowerlimit_Y = 0.45;
% upperlimit_Y = 0.53;

index_target_X = find(centers(:,1)>Xpixel*lowerlimit_X & centers(:,1)<Xpixel*upperlimit_X); % sort of handpicked targeted points...
index_target_Y = find(centers(:,2)>Ypixel*lowerlimit_Y & centers(:,2)<Ypixel*upperlimit_Y); 

index_target = intersect(index_target_X,index_target_Y);

radii_target = radii(index_target);
centers_target = centers(index_target,:); % [Xpixel Ypixel]
h_target = viscircles(centers_target,radii_target);

% position_error = [centers_target(1)/Xpixel centers_target(2)/Ypixel]

%% Draw lines between the targeted points
position = [0 0 0 0; centers_target(1,:) centers_target(4,:); centers_target(2,:) centers_target(3,:)]; % for some reason, the first row doesn't affect [x1 y1 x2 y2] which connects the two points (x1,y1) and (x2,y2)
image_lines = insertShape(image_filtered,'Line', position, 'Color', 'red','LineWidth',6);

figure(3); imshow(image_lines); title('Image with targeted lines'); 
h = viscircles(centers_target,radii_target);

%% Find rectangles on the image using 'regionprops'

bw = imbinarize(image); % Making the image as a binary

% find both black and white regions
stats = [regionprops(bw); regionprops(not(bw))]; % Find the objects on the image
CentroidDouble = cell2mat({stats.Centroid}.');

% Centroid limit
lowerlimitXpixel = 0.45;
upperlimitXpixel = 0.65;
lowerlimitYpixel = 0.45;
upperlimitYpixel = 0.65;

indexcentroid = find(CentroidDouble(:,1) > Xpixel*lowerlimitXpixel & CentroidDouble(:,1) < Xpixel*upperlimitXpixel & CentroidDouble(:,2) > Ypixel*lowerlimitYpixel & CentroidDouble(:,2) < Ypixel*upperlimitYpixel);

stats = stats(indexcentroid);

% Area limit
AreaDouble = cell2mat({stats.Area}.');

lowerlimitarea = 50000;
upperlimitarea = 400000;

indexarea = find(AreaDouble > lowerlimitarea & AreaDouble < upperlimitarea);
stats = stats(indexarea);

% show the image and draw the detected rectangles on it
figure(4); imshow(bw); hold on;

for i = 1:numel(stats)
    rectangle('Position',stats(i).BoundingBox,'Linewidth',3,'EdgeColor','y','LineStyle','-');
end


%% Collect the centered data

% 0504 camera-projector alignment measurement
Center_Test1 = [1550.37891508796,1024.81465684077];
Error_Test1 =[0.502065710844546,0.496518729089522];
percenterror_Test1 = (Center_Test1-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

Center_Test2 = [1546.66993855691,1046.05359163578];
Error_Test2 = [0.500864617408326,0.506808910676251];
percenterror_Test2 = (Center_Test2-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

Center_DMD = [1560.84359411161,1004.41406858377];
Error_DMD = [0.505454531771893,0.486634723151051];
percenterror_DMD = (Center_DMD-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

Center_Rect = [1560.83749552052,1003.85960177388];
Boundingbox_Rect = [1362.50000000000,889.500000000000,397,229]; % [x y Width Height]
% Xpixel_Rect = 
% Ypixel_Rect = 
percenterror_Rect = (Center_Rect-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

Centers_Flange_All = [944.947003986258,400.391219475034;2175.51511672542,406.487224684696;940.629598922905,1626.36327451950;2164.19667210200,1635.42500119377]
Center_Flange = sum(Centers_Flange_All)./4;
percenterror_Flange = (Center_Flange-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

% percenterror_Test1 = (Error_Test1-[0.5 0.5]).*100;
% percenterror_Test2 = (Error_Test2-[0.5 0.5]).*100;
% percenterror_DMD = (Error_DMD-[0.5 0.5]).*100;

% image = imread('Test1.tiff','tif');
figure(4); imshow(image); title('Original Image'); hold on;

plot(Center_Test1(1),Center_Test1(2),'r.');
plot(Center_Test2(1),Center_Test2(2),'r*');
plot(Center_DMD(1),Center_DMD(2),'g.');
plot(Center_Rect(1),Center_Rect(2),'g*');
plot(Center_Flange(1),Center_Flange(2),'c.');

xlabel('Xpixel (3088)')
ylabel('Ypixel (2064)')
legend('Position 1 (image)','Position 2','DMD circle','DMD rectangle','Flange plate');

% % Far distance
% center_TEST1_centercircle = [1559.23332014588,1100.00063956945] ;
% center_TEST1_bigcircle = [1560.51376473021,1120.82354756612]; 
% center_TEST1_rect = [1559.63832729212,1100.74869875258];
% 
% % Close distance
% center_TEST3_centercircle = [1554.12183847556,1105.36804894903];
% center_TEST3_rect = [1554.82470534931,1106.53393202558];
% 
% % LCPA1 center point from close distance 1 (11-inch)
% center_LCPA1_1 = [1535.95318335229,1119.68697763055];
% center_LCPA1_3 = [1535.69668863347,1118.67123674568];
% center_LCPA1_5 = [1535.93385411326,1119.51494830992];
% 
% % LCPA1 center point from close distance 2 (8.9-inch)
% center_LCPA1_2 = [1545.13931906187,1135.37225903406];
% center_LCPA1_4 = [1546.84754649983,1135.67557373299];
% center_LCPA1_6 = [1543.55845264212,1135.55767731425];
% 
% radii_LCAP1_1 = [18.9653379860763];
% 
% % Plot all the centers on the image
% figure(4); imshow(image); title('Original Image'); hold on;
% plot(center_TEST1_centercircle(1),center_TEST1_centercircle(2),'r.');
% plot(center_TEST1_bigcircle(1),center_TEST1_bigcircle(2),'r*');
% plot(center_TEST1_rect(1),center_TEST1_rect(2),'y+');
% 
% plot(center_TEST3_centercircle(1),center_TEST3_centercircle(2),'r.');
% plot(center_TEST3_rect(1),center_TEST3_rect(2),'y+');
% 
% plot(center_LCPA1_1(1),center_LCPA1_1(2),'c.'); 
% plot(center_LCPA1_3(1),center_LCPA1_3(2),'c.');
% plot(center_LCPA1_5(1),center_LCPA1_5(2),'c.');
% 
% plot(center_LCPA1_2(1),center_LCPA1_2(2),'c*'); 
% plot(center_LCPA1_4(1),center_LCPA1_4(2),'c*');
% plot(center_LCPA1_6(1),center_LCPA1_6(2),'c*');
% 
% xlabel('Xpixel (3088)')
% ylabel('Ypixel (2064)')
% legend('center circle (test pattern)','big circle (test pattern)','rectangle (test pattern)','','','Far distance (LCPA1)','','','Close distance (LCPA1)');

%% Find circles on the image using 'regionprops'

% bw = imbinarize(image);
% stats = regionprops('table',bw,'Centroid','MajorAxisLength','MinorAxisLength')
% 
% centers = stats.Centroid;
% diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
% radii = diameters/2;
% 
% hold on
% figure; imshow(bw); 
% viscircles(centers,radii);
% hold off

%% Further image conversions

% %% Erode
% SE = strel('disk',50,8);
% theImageEroded = imerode(image,SE);
% figure; imshow(theImageEroded); title('Eroded Image');

% %% Binarize image
% threshold = 25;
% [nRows,nCols] = size(testImage_1_);
% for ii = 1:nRows
%     for jj = 1:nCols
%         if (testImage_1_(ii,jj) > threshold)
%             binaryImage(ii,jj) = true;
%         else
%             binaryImage(ii,jj) = false;
%         end
%     end
% end
