% CheckSpectralReproducibility.
% 
% This is a simple program for checking spectral reproducibility over
% setting the trombone / bite-bar assembly positions.
%
% We observed the differences on spectral power even when setting the
% exactly the same position of bite-bar assembly and trombone positions.
% Therefore, this is to check how deviated the spectral results are by
% setting the same settings repeatedly.

% History:
%    12/22/22   smo   - Wrote it.

%% Initialize.
clear; close all;

%% Set test filename to save.
testFilename = 'Trombone';

%% Put image on the screen.
%
% Here we use plain white screen and channel settings are all turned on as
% its maximum.
WarmupScreen;

%% Open spectroradiometer.
OpenSpectroradiometer;

%% Measurement happens here.
% 
% Press any key before making a measurement. Every before making a key
% press, we will reset the Trombone setting (or the setting that you're
% interested). Basically the measurement will be done with the exactly the
% same settings, but we are interested in how reproducible it is.

nMeasures = 10;
for mm = 1:nMeasures
    % Key stroke to start measurement.
    disp('Press a key to start measurement!');
    pause;
    
    % Measure it. 
    spd(:,mm) = MeasureSpectroradiometer;
    
    % Show the progress.
    fprintf('Measurement Progress %d/%d \n', mm, nMeasures);    
end

%% Save the results.
SAVERESULTS = false;
if (SAVERESULTS)
    if (ispref('SpatioSpectralStimulator','SACCMaterials'))
        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCMaterials'),'CheckSpectralReproducibility');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilenameMeasure = sprintf('%s_%s',testFilename,dayTimestr);
        save(testFilenameMeasure,'spd');
    end
end

%% Calculate the scale factor.
%
% We set one spectrum as a reference to calculate the scale factor.
numSpdRef = 5;
numSpdTest = setdiff([1:1:size(spd,2)], numSpdRef);
spdRef = spd(:,numSpdRef);
spdTest = spd(:,numSpdTest);

% Calculate the scale factor here.
nSpdTest = length(numSpdTest);
for tt = 1:nSpdTest
    scaleFactor(tt) = spdTest(:,tt)\spdRef;
end

% Get mean, min, and max scale factors.
minScaleFactor = min(scaleFactor);
maxScaleFactor = max(scaleFactor);
meanScaleFactor = mean(scaleFactor);

%% Plot the results.
%
% Set the wavelength range.
S = [380 2 201];
wls = SToWls(S);

% Plot it.
figure; hold on;
plot(wls,spd(:,numSpdRef),'r-','linewidth',2);
plot(wls,spd(:,numSpdTest),'k-');
xlabel('Wavelength (nm)','fontsize',12);
ylabel('Spectral output (no unit)','fontsize',12);
xlim([380 780]);
title('Spectral Reproducibility Results','fontsize',12);
legend('Ref','Test');

% Add scale factor info to the plot.
main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
text(0.15, 0.9, sprintf('Mean scale factor = (%.2f)',meanScaleFactor), 'Parent', main,'fontsize',12)            
text(0.15, 0.8, sprintf('Min scale factor   = (%.2f)',minScaleFactor), 'Parent', main,'fontsize',12)
text(0.15, 0.75, sprintf('Max scale factor  = (%.2f)',maxScaleFactor), 'Parent', main,'fontsize',12)

%% Save the plot.
SAVETHEPLOT = true;
if (SAVETHEPLOT)
    testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCMaterials'),'CheckSpectralReproducibility');
    testFilenamePlot = fullfile(testFiledir, testFilename);
    testFileFormat = '.tiff';
    saveas(gcf,append(testFilenamePlot,testFileFormat));
    fprintf('\t Plot has been saved successfully! \n');
end

%% Close.
CloseSpectroradiometer;
CloseSceen;
