% SACC_FitCSFReadOnly
%
% This is to read the existing fitting results.
%
% See also:
%    SACC_FitCSFReadOnly.

% History:
%    3/28/23   smo    - Wrote it.
%    4/24/23   smo    - Now we save the CCSF figure and also read the
%                       results continuosly by key press.

%% Initialize.
clear; close all;

%% Set variables.
SaveFigure = true;
PlotAllSubjects = true;
WaitForKeyPress = true;

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
        
        % Check if each subject has CSF plots for all filters.
        nFiles = 0;
        for ff = 1:length(filterOptions)
            fileNameTemp = sprintf('%s_%s.%s',targetFilename,filterOptions{ff},imgFormat);
            for ll = 1:length(fileList)
                nFiles = nFiles + contains(fileNameTemp,fileList(ll).name);
            end
        end
        
        % Save if the subject has all available figures.
        if (nFiles == length(filterOptions))
            subjectAvailable{end+1} = whichSub;
        end
    end
end

%% Prompt which subject result to read.
%
% If we are not plotting all subjects, it will ask you which subject result
% to plot.
if (~PlotAllSubjects)
    while 1
        inputMessage = sprintf('\t Choose a subject to see the results: \n \t [%s] \n',strjoin(subjectAvailable));
        whichSub = input(inputMessage, 's');
        
        if ismember(whichSub, subjectAvailable)
            break
        end
        
        disp('Choose a number within the numbers displaying!');
    end
    
    fprintf('Subject (%s) results will be displaying...\n',whichSub);
end

%% Display the results.
if (PlotAllSubjects)
    nSubjects = length(subjectAvailable);
else
    nSubjects = 1;
end

for ss = 1:nSubjects
    % Choose subject to test if running for all.
    if (PlotAllSubjects)
        whichSub = subjectAvailable{ss};
    end
    
    % Get file location.
    testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
    testFiledir = fullfile(testFiledir,whichSub,'CSF');
    
    % Make a figure window.
    figure;
    set(gcf,'position',figurePosition);
    
    % Make a loop to load all results of the filters.
    for ff = 1:length(filterOptions)
        % Make a subplot for each filter.
        subplot(2,3,ff);
        
        % Get the filename to plot.
        testFilename = fullfile(testFiledir,...
            sprintf('CSF_%s_%s.%s',whichSub,filterOptions{ff},imgFormat));
        image = imread(testFilename);
        
        % Resize the image if you want.
        imgMagnifyIndex = 1;
        image = imresize(image,imgMagnifyIndex);
        [imgSizeHorz, imgSizeVert, ~] = size(image);
        
        % Display the image.
        % Crop image to remove the legend which is outside the plot.
        imageCropped = imcrop(image,[0 0 imgSizeVert-625 imgSizeHorz]);
        imshow(imageCropped);
        
        % Add title per each image.
        title(sprintf('Filter %s',filterOptions{ff}),'fontsize',20);
        
        % Add legend.
        if ff == length(filterOptions)
            subplot(2,3,ff+1);
            imageLegend = imcrop(image,[1040 90 525 150]);
            imshow(imageLegend);
        end
    end
    
    % Save the results.
    if (SaveFigure)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = getpref('SpatioSpectralStimulator','SACCAnalysis');
            testFiledir = fullfile(testFiledir,whichSub,'CSF');
            testFilename = fullfile(testFiledir,...
                sprintf('CSF_%s_AllFiltersSmooth.%s',whichSub,imgFormat));
            saveas(gcf, testFilename);
            fprintf('\t CSF plot has been successfully saved! \n');
        end
    end
    
    % Key press to draw next plot.
    if (PlotAllSubjects)
        % Show the progess.
        fprintf('\t Press a key to draw next plot! - (%d/%d) \n',ss,nSubjects);
        
        % Wait for the key press to plot next one.
        if (WaitForKeyPress)
            pause;
        end
        
        % Close the figure.
        close all;
    end
end
