% t_CSFGeneratorCheckTestImages
%
% This is to check test images quickly by eye.

% History:
%    10/25/22   smo    - Wrote it.

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
end

%% Plot images.
figure; clf;

for ss = 1:nSineFreqCyclesPerDeg
    sineFreqCyclesPerDegTemp = sineFreqCyclesPerDeg(ss);
  
    subplot(1,nSineFreqCyclesPerDeg,ss);
    imshow(images{ss});
    title(sprintf('%dcpd',sineFreqCyclesPerDegTemp));
end