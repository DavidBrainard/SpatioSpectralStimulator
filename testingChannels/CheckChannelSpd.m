% CheckChannelSpd
%
% This is to check channel power measuremtns.
%
% See Also:
%    MeasureChannelSpd.

% History:
%    03/29/22 dhb, smo  - Add in analysis.
%    03/31/22 smo       - Made the part of calculation of k as a function.
%    05/05/22 smo       - Added an option to load power meter data in a
%                         general file name.

%% Initialize.
clear; close all;

%% Set parameters.
nPrimaries = 3;
projectorModeNormal = false;
powerMeterWl = 550;
VERBOSE = true;

% Set this number to 0 if you are not willing to use Dataset 1,2,3. 
POWERDATASET = 0;

%% Load spectrum data here.
DEVICE = 'PR670';

% Make a string for file name.
switch projectorModeNormal
    case true
        projectorMode = 'NormalMode';
    case false
        projectorMode = 'SteadyOnMode';
end

% We will load different (black corrected or not) spd data over the power
% meter data set.
%
% Measurement date 3/28 is no black correction,
% 3/30 and 4/19 are black corrected.
switch POWERDATASET
    case 1
        olderDate = 2;
    case 2
        olderDate = 2;
    case 3
        olderDate = 0;
    otherwise
        olderDate = 0;
end

% Load the data here.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    testFilename = GetMostRecentFileName(testFiledir,sprintf('SpdData_%s_%s',DEVICE,projectorMode),'olderDate',olderDate);
    prData = load(testFilename);
else
    error('Cannot find data file');
end

% Cut the negative parts on the spectrum which caused by black correciton.
for pp = 1:nPrimaries
    prData.spdMeasured{pp} = max(prData.spdMeasured{pp},0);
end

% Set some params.
targetChannels = prData.targetChannels;
nTargetChannels = length(targetChannels);
S = prData.S;

%% Load powermeter data here.
curDir = pwd;
cd(testFiledir);

DEVICE = 'PowerMeter';
fileType = '.csv';

switch POWERDATASET
    case 1
        % DATASET 1.
        % Fixed wavelength sensitivity (550 nm).
        powerSingleNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalSingle');
        powerSingleSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnSingle');
        powerWhiteNormalWatt = xlsread('PowerMeterProcessedData.xlsx','NormalWhite');
        powerWhiteSteadyOnWatt = xlsread('PowerMeterProcessedData.xlsx','SteadyOnWhite');
        
        if (projectorModeNormal)
            powerMeterWatt = powerSingleNormalWatt;
            powerMeterWhiteWatt = powerWhiteNormalWatt;
        else
            powerMeterWatt = powerSingleSteadyOnWatt;
            powerMeterWhiteWatt = powerWhiteSteadyOnWatt;
        end
        
    case 2
        % DATASET 2.
        % Data with different wavelength.
        date = '0329';
        powerMeterWls = [448 476 404 552 592 620];
        dataRange = 'D17:D17';
        
        % Load sinlge peak data here.
        for pp = 1:nPrimaries
            for cc = 1:nTargetChannels
                targetChPeakWl = powerMeterWls(cc);
                targetCh = targetChannels(cc);
                fileName = append(DEVICE,'_',projectorMode,'_Primary',...
                    num2str(pp),'_Ch',num2str(targetCh),'_',num2str(targetChPeakWl),'nm_',date,fileType);
                readFile = readmatrix(fileName, 'Range', dataRange);
                powerMeterWatt(cc,pp) = readFile;
            end
        end
        
        % Load white data here.
        for cc = 1:nTargetChannels
            targetChPeakWl = powerMeterWls(cc);
            fileName = append(DEVICE,'_',projectorMode,'_White_',num2str(targetChPeakWl),'nm_',date,fileType);
            readFile = readmatrix(fileName, 'Range', dataRange);
            powerMeterWhiteWatt(cc,:) = readFile;
        end
        
    case 3
        % DATASET 3.
        % Black corrected.
        % Fixed wavelength sensitivity (550 nm).
        date = '0330';
        fileName = append(DEVICE,'_',projectorMode,'_Singles_',num2str(powerMeterWl),'nm_',date,fileType);
        readFile = readmatrix(fileName);
        powerMeterAllWatt = readFile;
        powerMeterWhiteWatt = powerMeterAllWatt(1,:);
        powerMeterWatt = powerMeterAllWatt(2:end,:);
        
    otherwise
        % Load the power meter data in general file name.
        nMeasurements = 19;
        dataRange = append('D16:D',num2str(16+nMeasurements-1));
        date = '0506';
        fileName = append(DEVICE,'_',projectorMode,'_',...
            num2str(powerMeterWl),'nm_',num2str(date),fileType);
        readFile = readmatrix(fileName, 'Range', dataRange);
        powerMeterWhiteWatt = readFile(1,:);
        powerMeterWatt = readFile(2:end,:);
        
end

% Match the power meter array size.
powerMeterWatt = reshape(powerMeterWatt,nTargetChannels,nPrimaries);

%% Plot measured spectra.
if (VERBOSE)
    % Single peak spectrum.
    for pp = 1:nPrimaries
        figure; clf;
        plot(SToWls(S),prData.spdMeasured{pp});
        title(append('Screen Primary: ',num2str(pp),' ',projectorMode),'FontSize',15);
        xlabel('Wavelength (nm)','FontSize',15);
        ylabel('Relative Spectral Power','FontSize',15);
        
        % Make strings for graph legend.
        for cc = 1:nTargetChannels
            targetChannel = prData.targetChannels(cc);
            legendSinglePeaks{cc} = append('Ch ',num2str(targetChannel));
        end
        legend(legendSinglePeaks);
    end
    
    % White.
    figure; clf;
    plot(SToWls(S),prData.spdMeasuredWhite);
    title(append('White ',projectorMode),'FontSize',15);
    xlabel('Wavelength (nm)','FontSize',15);
    ylabel('Relative Spectral Power','FontSize',15);
end

%% Find scale factors for each measurement
%
% Sinlge peaks.
for pp = 1:nPrimaries
    for cc = 1:nTargetChannels
        if (POWERDATASET == 2)
            powerMeterWl = powerMeterWls(cc);
        end
        k(cc,pp) = SpdToPower(prData.spdMeasured{pp}(:,cc), powerMeterWatt(cc,pp), 'targetWl', powerMeterWl)';
    end
end

% White.
nWhites = length(powerMeterWhiteWatt);
for ww = 1:nWhites
    if (POWERDATASET == 2)
        powerMeterWl = powerMeterWls(ww);
    end
    kWhite(ww) = SpdToPower(prData.spdMeasuredWhite, powerMeterWhiteWatt(ww), 'targetWl', powerMeterWl);
end

% Plot it.
if (VERBOSE)
    figure; clf; hold on;
    targetChannelsWhite = targetChannels(1:length(kWhite));
    plot(targetChannelsWhite,kWhite,'bo','MarkerSize',12,'MarkerFaceColor','b');
    plot(targetChannels,k,'ro','MarkerSize',12,'MarkerFaceColor','r');
    xlabel('Target Channel','FontSize',15);
    ylabel('Conversion Coefficient k','FontSize',15);
    xticks(targetChannels);
    ylim([0 0.03]);
    legend('White','Single peak','FontSize',13);
    title(append('DataSet ',num2str(POWERDATASET),' ',projectorMode),'FontSize',15);
end
