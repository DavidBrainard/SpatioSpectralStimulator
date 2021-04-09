
%% Initialize
clear; close all;

%% Load in an image
theImage = imread('~/Desktop/testImage.tiff','tif');
figure; imshow(theImage); title('Original Image');

%% Find size of circles we want
%d = imdistline;
theImageFilt = imfilter(theImage,ones(5)/25);
[centers, radii] = imfindcircles(theImageFilt,[110 140]/2,'ObjectPolarity','dark','Sensitivity',0.92,'EdgeThreshold',0.05);
h = viscircles(centers,radii);

%% Erode
SE = strel('disk',50,8);
theImageEroded = imerode(theImage,SE);
figure; imshow(theImageEroded); title('Eroded Image');

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
