% MakeContrastImage.
%
% This code makes vertical black and white stripe contrast image. 

% History:
%    08/14/23   smo   Wrote it.

%% Define image dimensions
imageWidth = 1920;
imageHeight = 1080;

% Create a blank image with white background
image = ones(imageHeight, imageWidth);

% Define the number of stripes and total stripe width
%
% Set numStrips = 50 (18 cpd) / 34 (12 cpd) / 25 (9 cpd) 17 (6 cpd) / 9 (3
% cpd). This is the number when measuring the printed pattern outside the
% SACCSFA room, on the hallway of 4th floor at Goddard.
numStripes = 9;
totalStripeWidth = imageWidth / numStripes;

%% Generate the alternating black and white stripes
for i = 1:numStripes
    startColumn = round((i - 1) * totalStripeWidth) + 1;
    endColumn = round(i * totalStripeWidth);
    
    % Set the pixels within the black stripe to black
    if mod(i, 2) == 1
        image(:, startColumn:endColumn) = 0;
    end
end

%% Display the image
figure;
imshow(image);

%% Save the image.
SAVEIMAGE = true;
if (SAVEIMAGE)
    imwrite(image,'image.tiff');
end