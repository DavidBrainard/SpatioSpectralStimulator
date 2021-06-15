%% Initialize
clear; close all;

%% Load an image and find targeted circle

% for i = 1:15;
% BaseName='0512_position';
% FileName=[BaseName,num2str(i)]

FileName=['0609 blkandwhite'];

% image = imread('Repeatability_position1_1.tiff','tif');
image = imread(FileName,'tif');
[Ypixel Xpixel] = size(image);

figure(1); imshow(image); title('Original Image');

% Find size of circles we want
d = imdistline; % measure the radius of the target circle image

image_filtered = imfilter(image,ones(5)/25);

[centers, radii] = imfindcircles(image_filtered,[64 80]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
% 60 80

h = viscircles(centers,radii);

% Find targeted circles 
figure(2); imshow(image_filtered); title('Image with targeted points'); hold on;

% % Centered circle on the test pattern
lowerlimit_X = 0.1;
upperlimit_X = 0.5;
lowerlimit_Y = 0.1;
upperlimit_Y = 0.99;

% % LCPA1 centered point 
% lowerlimit_X = 0.48;
% upperlimit_X = 0.51;
% lowerlimit_Y = 0.48;
% upperlimit_Y = 0.52;

index_target_X = find(centers(:,1)>Xpixel*lowerlimit_X & centers(:,1)<Xpixel*upperlimit_X); % sort of handpicked targeted points...
index_target_Y = find(centers(:,2)>Ypixel*lowerlimit_Y & centers(:,2)<Ypixel*upperlimit_Y); 

index_target = intersect(index_target_X,index_target_Y);

radii_target = radii(index_target);
centers_target = centers(index_target,:); % [Xpixel Ypixel]
h_target = viscircles(centers_target,radii_target);

% centers_target_collect(i,:) = centers_target
% position_error_collect(i,:) = [centers_target(1)/Xpixel centers_target(2)/Ypixel]

% end
%% Draw lines between the targeted points
% position = [0 0 0 0; centers_target(1,:) centers_target(4,:); centers_target(2,:) centers_target(3,:)]; % for some reason, the first row doesn't affect [x1 y1 x2 y2] which connects the two points (x1,y1) and (x2,y2)
% image_lines = insertShape(image_filtered,'Line', position, 'Color', 'red','LineWidth',6);
% 
% figure(3); imshow(image_lines); title('Image with targeted lines'); 
% h = viscircles(centers_target,radii_target);

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

center_target_rect = stats.Centroid;
center_boundingbox = stats.BoundingBox;

%% Collect the centered data

%% 0518 Alignment offset data
% Center data
centers_pin = [697.391479038494,1978.69204932857; 2495.92033516091,221.416817516480];
raddi_pin = [35.9888791200350; 37.8473222670018];

center_pin_left = centers_pin(1,:);
center_pin_right = centers_pin(2,:);
raddi_pin_left = raddi_pin(1,:);
raddi_pin_right = raddi_pin(2,:);

center_DMD = [1546.66429119242,1025.86867106685];
center_rect = [1547.02309342881,1025.81268004695];

center_boundingbox = [1351.50000000000,915.500000000000,390,221];
size_DMD = center_boundingbox(3:4);
size_DMD_diagonal = sqrt(size_DMD(1)^2+size_DMD(2)^2);

% Offset data (in pixel unit)
offset_left = center_pin_left-center_DMD;
offset_left_diagonal = sqrt(offset_left(1)^2+offset_left(2)^2);
offset_right = center_pin_right-center_DMD;
offset_right_diagonal = sqrt(offset_right(1)^2+offset_right(2)^2);

% Flange plate data
centers_flange=[933.513373804276,416.481730715954;931.181518664654,1645.33208597354;2162.50249355838,423.012431801115;2149.48977613578,1649.68807857790]

length_flange1 = pdist([centers_flange(1,:);centers_flange(2,:)],'euclidean')
length_flange2 = pdist([centers_flange(2,:);centers_flange(4,:)],'euclidean')
length_flange3 = pdist([centers_flange(3,:);centers_flange(4,:)],'euclidean')
length_flange4 = pdist([centers_flange(1,:);centers_flange(3,:)],'euclidean')

length_flange5 = pdist([centers_flange(1,:);centers_flange(4,:)],'euclidean') % diagonal
length_flange6 = pdist([centers_flange(2,:);centers_flange(3,:)],'euclidean') % diagonal 

lengths_flange = [length_flange1; length_flange2; length_flange3; length_flange4]
lengths_flange_diagonal = [length_flange5; length_flange6]

% Center pins location displacement
centers_threepins=[650.486120674839,127.064398769850; 658.867060372274,1859.45761796232; 2415.94558949847,112.070167735241];

offset_threepins = centers_threepins - center_DMD;
offset_threepins_diagonal = sqrt(offset_threepins(:,1).^2+offset_threepins(:,2).^2)

% Pixeltomm
pixeltomm_DMD = 0.0535398305;
pixeltomm_Flange = sqrt(60^2+60^2)/mean(lengths_flange_diagonal);
pixeltomm_calipers = 128.12 /(raddi_pin(1)+raddi_pin(2)+pdist([centers_pin(1,:);centers_pin(2,:)],'euclidean'));% 128.12 mm is the outer distance between two pins
pixeltomm_Freddie = 125.35 / pdist([centers_pin(1,:);centers_pin(2,:)],'euclidean'); % Pins' center to center distance

pixeltomm = pixeltomm_Freddie; % Set the conversion factor here

offset_left_mm = offset_left.*pixeltomm;
offset_left_diagonal_mm = offset_left_diagonal.*pixeltomm;
offset_right_mm = offset_right.*pixeltomm;
offset_right_diagonal_mm = offset_right_diagonal.*pixeltomm;

offset_threepins_mm = offset_threepins.*pixeltomm; % [top left; bottom left; top right]
offset_threepins_diagoanl_mm = offset_threepins_diagonal.*pixeltomm;

% figure; imshow(image_filtered); title('Image with targeted points'); hold on;
% h_target = viscircles(centers_pin,raddi_pin);

%% 0512 Camera center alignment data (5-inch apart) after fixing the camera
% plate

% image = imread('0512_position1.tiff','tif');
% [Ypixel Xpixel] = size(image);
% 
center_position1=[1546.03236177467,1031.48147239997];
percenterror_position1 = (center_position1-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

center_position2=[1543.88214884096,1031.71576750248];
percenterror_position2 = (center_position2-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

center_DMD=[1547.31722148188,1027.10083644431];
percenterror_DMD = (center_DMD-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

center_rect=[1548.20250417750,1027.78683160045];
percenterror_rect = (center_rect-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

centers_flange=[933.595599705120,417.797163118484;2163.59504090741,423.398456555532;2151.53488595093,1652.50448137099;928.789883108596,1647.79719001034];
center_flange=mean(centers_flange);
percenterror_flange = (center_flange-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

errormm_position2 = (center_position2 - center_position1)*PixelTomm_DMD;
errormm_DMD = (center_DMD - center_position1)*PixelTomm_DMD;
errormm_rect = (center_rect - center_position1)*PixelTomm_DMD;
errormm_flange = (center_flange - center_position1)*PixelTomm_DMD;

centers_threepins=[650.486120674839,127.064398769850; 658.867060372274,1859.45761796232; 2415.94558949847,112.070167735241];
% bottom left from the small hole (offset from the center) = [666.427387394665,1848.24547192977]

% figure(4); imshow(image); title('Original Image'); hold on;
% plot(center_position1(1), center_position1(2),'r.');
% plot(center_position2(1), center_position2(2),'r*');
% plot(center_DMD(1), center_DMD(2),'g.');
% plot(center_rect(1), center_rect(2),'g*');
% plot(center_flange(1), center_flange(2),'y+');
% 
% legend('Position1(image)','Position2','DMD circle','DMD rect','Flange plate');
% xlabel('Xpixel (3088)')
% ylabel('Ypixel (2064)')


%% 0506 Camera repeatability test 
% load('Position1_Center_Pixel.mat')
% load('Position2_Center_Pixel.mat')
% 
% centers_target_collect_1 = centers_target_collect_1(6:15,:);
% centers_target_collect_2 = centers_target_collect_2(6:15,:);
% 
% mean_position1 = mean(centers_target_collect_1);
% mean_position2 = mean(centers_target_collect_2);
% std_position1 = std(centers_target_collect_1);
% std_position2 = std(centers_target_collect_2);
% 
% max_position1_X = max(centers_target_collect_1(:,1))-min(centers_target_collect_1(:,1));
% max_position1_Y = max(centers_target_collect_1(:,2))-min(centers_target_collect_1(:,2));
% max_position1 = [max_position1_X max_position1_Y];
% 
% max_position2_X = max(centers_target_collect_2(:,1))-min(centers_target_collect_2(:,1));
% max_position2_Y = max(centers_target_collect_2(:,2))-min(centers_target_collect_2(:,2));
% max_position2 = [max_position2_X max_position2_Y];


%% 0504 camera-projector alignment measurement
% Center_Test1 = [1550.37891508796,1024.81465684077];
% Error_Test1 =[0.502065710844546,0.496518729089522];
% percenterror_Test1 = (Center_Test1-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;
% 
% Center_Test2 = [1546.66993855691,1046.05359163578];
% Error_Test2 = [0.500864617408326,0.506808910676251];
% percenterror_Test2 = (Center_Test2-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;
% 
% Center_DMD = [1560.84359411161,1004.41406858377];
% Error_DMD = [0.505454531771893,0.486634723151051];
% percenterror_DMD = (Center_DMD-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;
% 
% Center_Rect = [1560.83749552052,1003.85960177388];
% Boundingbox_Rect = [1362.50000000000,889.500000000000,397,229]; % [x y Width Height]
% % Xpixel_Rect = 
% % Ypixel_Rect = 
% percenterror_Rect = (Center_Rect-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;
% 
% Centers_Flange_All = [944.947003986258,400.391219475034;2175.51511672542,406.487224684696;940.629598922905,1626.36327451950;2164.19667210200,1635.42500119377]
% Center_Flange = sum(Centers_Flange_All)./4;
% percenterror_Flange = (Center_Flange-0.5*[Xpixel Ypixel])./[Xpixel Ypixel]*100;

% image = imread('Test1.tiff','tif');
% figure(4); imshow(image); title('Original Image'); hold on;
% plot(centers_target_collect_1(:,1), centers_target_collect_1(:,2),'r.');
% plot(centers_target_collect_2(:,1), centers_target_collect_2(:,2),'g.');
% legend('Position1(N=10)','Position2(N=10)');
% xlabel('Xpixel (3088)')
% ylabel('Ypixel (2064)')

% 
% plot(Center_Test1(1),Center_Test1(2),'r.');
% plot(Center_Test2(1),Center_Test2(2),'r*');
% plot(Center_DMD(1),Center_DMD(2),'g.');
% plot(Center_Rect(1),Center_Rect(2),'g*');
% plot(Center_Flange(1),Center_Flange(2),'c.');
% 
% xlabel('Xpixel (3088)')
% ylabel('Ypixel (2064)')
% legend('Position 1 (image)','Position 2','DMD circle','DMD rectangle','Flange plate');

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
