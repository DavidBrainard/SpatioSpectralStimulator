% SACC_analyzeMonitor.
%
% This is for quick check of the calibration results.

%% Initialize.
clear; close all;

%% Set dir.
dir = cd;

fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Calibration';
cd(fileDir);

%% Load data.
dataSACC = load('SACC.mat').cals;
dataSACCPrimary1 = load('SACCPrimary1.mat').cals;
dataSACCPrimary2 = load('SACCPrimary2.mat').cals;
dataSACCPrimary3 = load('SACCPrimary3.mat').cals;

% Get the most recent data.
nPrimaries = 3;
numMostRecentDataSACC = length(dataSACC);
for pp = 1:nPrimaries
    numMostRecentDataSACCPrimary(pp) = length(eval(append('dataSACCPrimary',num2str(pp))));
end

nData = 2;
for dd = 1:nData
    numTargetSACC(dd) = numMostRecentDataSACC - dd + 1;
    numTargetSACCPrimary(:,dd) = numMostRecentDataSACCPrimary - dd + 1;
    
    % Spd.
    spdSACC{dd} = dataSACC{numTargetSACC(dd)}.processedData.P_device;
    spdPrimary1{dd} = dataSACCPrimary1{numTargetSACCPrimary(1,dd)}.processedData.P_device;
    spdPrimary2{dd} = dataSACCPrimary2{numTargetSACCPrimary(2,dd)}.processedData.P_device;
    spdPrimary3{dd} = dataSACCPrimary3{numTargetSACCPrimary(3,dd)}.processedData.P_device;
    
    % Date of measurement.
    dateSACC{dd} = dataSACC{numTargetSACC(dd)}.describe.date;
    dateSACCPrimary1{dd} = dataSACCPrimary1{numTargetSACCPrimary(1,dd)}.describe.date;
    dateSACCPrimary2{dd} = dataSACCPrimary2{numTargetSACCPrimary(2,dd)}.describe.date;
    dateSACCPrimary3{dd} = dataSACCPrimary3{numTargetSACCPrimary(3,dd)}.describe.date;
end

S = [380 2 201];
wls = SToWls(S);

%% Plot it.
figure; hold on;
lineColor = {'r--','k-'};
lineWidth = 1;
fontSizeTitle = 15;
fontSizeAxis = 12;

% SACC
subplot(4,1,1); hold on;
for dd = 1:nData
    plot(wls,spdSACC{dd},lineColor{dd},'linewidth',lineWidth);
end
xlabel('Wavelength (nm)', 'fontsize', fontSizeAxis);
ylabel('Spectral power','fontsize', fontSizeAxis);
title('SACC','fontsize',fontSizeTitle);
legend(dateSACC);

% SACC primaries
for pp = 1:nPrimaries
    subplot(4,1,pp+1); hold on;
    xlabel('Wavelength (nm)','fontsize', fontSizeAxis);
    ylabel('Spectral power','fontsize', fontSizeAxis);
    title(append('SACC Primary ',num2str(pp)),'fontsize',fontSizeTitle);
    
    for dd = 1:nData
        spdTemp = eval(append('spdPrimary',num2str(pp)));
        plot(wls,spdTemp{dd},lineColor{dd},'linewidth',lineWidth);
    end
    legend(eval(append('dateSACCPrimary',num2str(pp))));
end

%% Back to current dir.
cd(dir);
