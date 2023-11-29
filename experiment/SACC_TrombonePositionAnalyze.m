% SACC_TrombonePositionAnalyze
%
% This simply shows the distribution of the Trombone positions of all
% subjects in the SACC experiment.
%
% See also:
%    SACC_calibrateTrombone.

% History:
%    11/29/23   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Read out the experimental log table.
testFiledir = getpref('SpatioSpectralStimulator','SACCMaterials');
testFiledir = fullfile(testFiledir,'Experiment','ExperimentalLog');
testFilename = 'SACC_Experimental_Log.xlsx';
testFileToRead = fullfile(testFiledir,testFilename);

% Get available sheet names of the file.
[~, sheetNames] = xlsfinfo(testFileToRead);

% Get the subject options that completed the study.
nSheets = length(sheetNames);
subjectOptions = {};
trombone = [];
for ss = 1:nSheets
    tableTemp = readtable(testFileToRead,'sheet',sheetNames{ss});
    numVisitComplete = 7;
    if ismember('VisitNo_',tableTemp.Properties.VariableNames)
        if max(tableTemp.VisitNo_ == numVisitComplete)
            subjectOptions{end+1} = tableTemp.SubjectID(1);
            trombone(end+1) = tableTemp.Trombone(1);
        end
    end
end

%% Convert trombone position to refraction.
%
% Polyfit values ('p') from our measurement data. This is the fitting
% parameters when fitting Trombone position (mm) over Testing lens diopters
% within the range from -6 to +6. FYI, if a subject has refraction of -2,
% we need the lens diopter of +2 to correct it. So, the sign of the testing
% lens is is inverse of the subject's refraction.
p = [6.2225  150.9583];

% Get the trombone position of emmentropic.
lensDiopterEmmentropic = 0;
tromboneEmmentropic = polyval(p,lensDiopterEmmentropic);

% Convert trombone position into the values of lens diopters and refraction.
lensDiopter = (trombone - p(2))./p(1);
refraction = -lensDiopter;

%% Plot the histogram of the trombone position.
figure; hold on;
nSubjects = length(subjectOptions);

% Histogram.
histData = histogram(trombone,'FaceColor','y');
xlabel('Trombone position (mm)','fontsize',15);
ylabel('Population','fontsize',15);
title(sprintf('Trombone position distribution (N = %d)',nSubjects),'fontsize',15);

% Add reference line of emmentropic.
plot([tromboneEmmentropic tromboneEmmentropic],[0 max(histData.Values)],...
    ':','color','k','linewidth',3);

% Add legend.
legend('All subjects','Emmentropic','fontsize',12);
