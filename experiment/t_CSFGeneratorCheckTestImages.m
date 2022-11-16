% t_CSFGeneratorCheckTestImages
%
% This is to check test images quickly by eye.

% History:
%    10/25/22   smo    - Wrote it.

%% Initialize.
clear; close all;
 
%% Load images.
sineFreqCyclesPerDeg = [3,6,9,12,18];
olderDate = 0;
nSineFreqCyclesPerDeg = numel(sineFreqCyclesPerDeg);

for ss = 1:nSineFreqCyclesPerDeg
    sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
    
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
        testFilename = GetMostRecentFileName(testFiledir,...
            sprintf('RunExpData_%d_cpd',sineFreqCyclesPerDegTemp),'olderDate',olderDate);
        theData{ss} = load(testFilename);
        images{ss} = theData{ss}.sceneParamsStruct.predefinedRGBImages{31};
    else
        error('Cannot find data file');
    end
    [filepath,filename,ext] = fileparts(testFilename);
    fprintf('\t *Image data has been successfully loaded! (%d/%d) \n',ss,nSineFreqCyclesPerDeg);
    fprintf('\t *Loaded file name = (%s) \n',filename);
end

%% Plot images.
for ss = 1:nSineFreqCyclesPerDeg
    sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
    figure; clf;
    imshow(images{ss});
    title(sprintf('%dcpd',sineFreqCyclesPerDegTemp),'fontsize',15);
end

%% Load images (high contrast).
sineFreqCyclesPerDeg = [18];
olderDate = 0;
nSineFreqCyclesPerDeg = numel(sineFreqCyclesPerDeg);

for ss = 1:nSineFreqCyclesPerDeg
    sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
    
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCData'),'TestImages');
        testFilename = GetMostRecentFileName(testFiledir,...
            sprintf('RunExpData_high_%d_cpd',sineFreqCyclesPerDegTemp),'olderDate',olderDate);
        theData{ss} = load(testFilename);
        images{ss} = theData{ss}.sceneParamsStruct.predefinedRGBImages{31};
    else
        error('Cannot find data file');
    end
    [filepath,filename,ext] = fileparts(testFilename);
    fprintf('\t *Image data has been successfully loaded! (%d/%d) \n',ss,nSineFreqCyclesPerDeg);
    fprintf('\t *Loaded file name = (%s) \n',filename);
end

%% Plot images.
for ss = 1:nSineFreqCyclesPerDeg
    sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
    figure; clf;
    imshow(images{ss});
    title(sprintf('%dcpd',sineFreqCyclesPerDegTemp),'fontsize',15);
end
