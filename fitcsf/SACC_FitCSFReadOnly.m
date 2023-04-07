% SACC_FitCSFReadOnly
%
% This is to read the existing fitting results.
%
% See also:
%    SACC_FitCSFReadOnly.

% History:
%    3/28/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set variables.
%
% If oneFigure set to true, we will plot all results on one figure.
% Otherwise, we will make a separate plot per each filter.
oneFigure = true;

fittingMode = 'crossValBootAcrossFmincon';
filterOptions = {'A', 'B', 'C', 'D', 'E'};
imgFormat = 'tiff';
figureSize = 2000;
figurePosition = [200 300 figureSize figureSize];

%% Check the subjects with all data.
if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFileList = dir(testFiledir);
    
    % Find available data of subject and spatial frequency.
    for tt = 1:length(testFileList)
        testFilenameList{tt}  = testFileList(tt).name;
    end
    idxSubjectName = find(str2double(testFilenameList)>0);
    subjectNameOptions = testFilenameList(idxSubjectName);
    
    % Get the subjects having all figures available.
    subjectAvailable = {};
    for ss = 1:length(subjectNameOptions)
        % Get the subject number.
        whichSub = subjectNameOptions{ss};
        
        % Set the directory of the figures saved per each subject.
        testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
        testFiledir = fullfile(testFiledir,whichSub,'CSF');
        targetFilename = append(sprintf('CSF_%s',whichSub));
        fileList = dir(append(fullfile(testFiledir,targetFilename),'*'));
        
        % Check the number of available figures per subject.
        nFiles = 0;
        for ff = 1:length(fileList)
            fileNameTemp = string(fileList(ff).name);
            % THIS PART SHOULD BE MODIFIED TO READ FILES PROPERLY
            nFiles = nFiles;
        end
        
        % Save if the subject has all available figures.
        if (nFiles == 5)
            subjectAvailable{end+1} = whichSub;
        end
    end
end

%% Prompt which subject result to read.
while 1
    inputMessage = sprintf('\t Choose a subject to see the results: \n \t [%s] \n',strjoin(subjectAvailable));
    whichSub = input(inputMessage, 's');
    
    if ismember(whichSub, subjectAvailable)
        break
    end
    
    disp('Choose a number within the numbers displaying!');
end

fprintf('Subject (%s) results will be displaying...\n',whichSub);

%% Display the results.
testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
testFiledir = fullfile(testFiledir,whichSub,'CSF');

% If we make only one plot to have all results.
if (oneFigure)
    figure;
    set(gcf,'position',figurePosition);
end

% Make a loop to load all results of the filters.
for ff = 1:length(filterOptions)
    % Draw everything on one figure or separate.
    if (oneFigure)
        subplot(2,3,ff);
    else
        figure;
        set(gcf,'position',figurePosition);
    end
    
    % Get the filename to plot.
    testFilename = fullfile(testFiledir,...
        sprintf('CSF_%s_%s.%s',whichSub,filterOptions{ff},imgFormat));
    image = imread(testFilename);
    
    % Resize the image if you want.
    imgMagnifyIndex = 1;
    image = imresize(image,imgMagnifyIndex);
    
    % Display the image.
    imshow(image);
    
    % Add title per each image.
    if (oneFigure)
        title(sprintf('Filter %s',filterOptions{ff}),'fontsize',20);
    end
end
