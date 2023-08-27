% SACC_ScreenStabilityCheck
%
% This tests screen stability for SACC project. It basically measures
% the plain screen repeatedely and
%
% History:
%    11/24/21 smo   Pulled out this part from the old code. It has been
%                     cleaned up using our SACC measurement functions.
%    01/07/22 smo   Added the option skipping the measurement.
%    03/24/22 smo   Added an option using nullImage in the SACC experiment
%                   to measure it over the time.

%% Initialize.
clear; close all;

%% Set parameters here.
%
% This code basically measures the spectrum of the plain screen
% automatically in a specific time interval.
%
% Set measurement time length and interval. The time is set in minute unit,
% which will be converted into second unit later in this code.
totalMeasureTimeMin = 120;
timeDelayBeforeEachMeasurementMin = 0.1;
timeDelayBeforeMeasurementSec = timeDelayBeforeEachMeasurementMin * 60;

% Measurement range.
S = [380 2 201];

% Verbose.
VERBOSE = true;
MEASURE = false;

%% Make a bit of time delay before the measurement starts.
% You can go out before it's done.
if (MEASURE)
    timeDelayGoOutSec = 3;
    fprintf('You have %d seconds to go out the room!',timeDelayGoOutSec);
    WaitSecs(timeDelayGoOutSec);
    
    %% Make screen and spectroradiometer ready.
    %
    % Open the plain screen. Simply set it as white here. It won't change
    % during the whole measurements.
    screenSettings = [0 0 0];
    [window windowRect] = OpenPlainScreen(screenSettings,'projectorMode',true,'verbose',VERBOSE);
    
    load('nullImage.mat');
    SetScreenImage(nullImage,window,windowRect);
    
    % Set channel settings.
    %     nChannels = 16;
    %     nPrimaries = 3;
    %     channelSettings = ones(nChannels,nPrimaries);
    load('screenPrimarySettings.mat');
    channelSettings = screenPrimarySettings;
    SetChannelSettings(channelSettings);
    
    % Connect to the spectroradiometer. We will use PR670 here.
    OpenSpectroradiometer;
    
    %% Measurements.
    %
    % The number of measurements are calculated based on the total measurement
    % time and time delay interval which were set from the above.
    %
    % +1 is the measurement at the cold state which is right after the screen turned on.
    nMeasurments = round(totalMeasureTimeMin / timeDelayBeforeEachMeasurementMin) + 1;
    allSpdMeasured = zeros(S(3),nMeasurments);
    
    % Measure it.
    for i=1:nMeasurments
        
        % Measurement happens here.
        allSpdMeasured(:,i) = MeasureSpectroradiometer;
        
        % Make a delay before next measurement.
        WaitSecs(timeDelayBeforeMeasurementSec);
    end
else
    % Load the data if the measurement is skipped.
    olderDate = 1;
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('stabilityCheck'),'olderDate',olderDate);
        load(testFilename);
    end
end

%% XYZ calculations
%
% Match the wavelength range.
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos, 683*T_xyzJuddVos, S);

% Calculate XYZ values.
XYZ = T_xyz * allSpdMeasured;
xyY = XYZToxyY(XYZ);

% Calculate color gamut too.
colorGamut = XYZToxyY(T_xyz);
colorGamut(:,end+1) = colorGamut(:,1);

%% Plot the data.
if (VERBOSE)
    labelFontSize = 25;
    titleFontSize = 25;
    
    % Spds.
    figure; clf; hold on;
    plot(SToWls(S),allSpdMeasured);
    xlabel('Wavelenth (nm)','fontsize',labelFontSize);
    ylabel('Spectral power distribution','fontsize',labelFontSize);
    
    % Luminance.
    if (~exist('nMeasurements'))
        nMeasurments = size(allSpdMeasured,2);
    end
    figure; clf;
    measurementTime = linspace(0, totalMeasureTimeMin, nMeasurments);
    plot(measurementTime, XYZ(3,:),'b.','markersize',10);
    xlabel('Time (min)', 'fontsize',labelFontSize);
    ylabel('Luminance (arb unit)', 'fontsize',labelFontSize);
    ylim([max(XYZ(3,:))*0.8 max(XYZ(3,:))*1.2]);
%     title('Luminance','fontsize',titleFontSize);
    
    % xy coordiantes.
    figure; hold on;
    plot(xyY(1,:), xyY(2,:), 'b.','markersize',10);
    plot(colorGamut(1,:),colorGamut(2,:),'k-');
    xlabel('CIE x','fontsize',labelFontSize);
    ylabel('CIE y','fontsize',labelFontSize);
    legend('Measurements','Spectral locus');
%     title('CIE xy chromaticity','fontsize',titleFontSize);
    
    % CIE x over the time.
    figure;
    plot(measurementTime,xyY(1,:),'b.','markersize',10);
    xlabel('Time (min)','fontsize',labelFontSize);
    ylabel('CIE x','fontsize',labelFontSize);
%     title('CIE x','fontsize',titleFontSize);
    
    % CIE y over the time.
    figure;
    plot(measurementTime,xyY(2,:),'b.','markersize',10);
    xlabel('Time (min)','fontsize',labelFontSize);
    ylabel('CIE y','fontsize',labelFontSize);
%     title('CIE y','fontsize',titleFontSize);
end

%% Save the data.
%
% Close the screen and projector.
if (MEASURE)
    CloseScreen;
    CloseSpectroradiometer;
    
    % Save data with the name containing dayTimestr.
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,sprintf('stabilityCheck_%s',dayTimestr));
        save(testFilename,'allSpdMeasured','XYZ','xyY','colorGamut','nMeasurments', ...
            'totalMeasureTimeMin','timeDelayBeforeEachMeasurementMin','S');
    end
end
