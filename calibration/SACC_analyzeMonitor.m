% SACC_analyzeMonitor.
%
% This is for quick check of the calibration results.

% History:
%    10/13/22   smo    - Wrote it.

%% Initialize.
clear; close all;

%% Set dir.
%
% Save the current dir. We will set it back to current dir after.
dir = cd;

% Get the dir that stors the measurement files.
fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/Calibration';
cd(fileDir);

%% Load calibration data.
dataSACC = load('SACC.mat').cals;
dataSACCPrimary1 = load('SACCPrimary1.mat').cals;
dataSACCPrimary2 = load('SACCPrimary2.mat').cals;
dataSACCPrimary3 = load('SACCPrimary3.mat').cals;

% Get a number for each data which contains the most recent data.
%
% SACC.
olderDate = 0;
numMostRecentDataSACC = length(dataSACC) - olderDate;

% SACC primaries.
nPrimaries = 3;
for pp = 1:nPrimaries
    numMostRecentDataSACCPrimary(pp) = length(eval(append('dataSACCPrimary',num2str(pp)))) - olderDate;
end

% Extract spd and the date of measurement from each dataset. You can choose
% multiple from each dataset if you want. It reads the data from the most
% recent one, so if you set nData = 2, it will read the most recent one and
% the second to the most recent.
nData = 3;
for dd = 1:nData
    numTargetDataSACC(dd) = numMostRecentDataSACC - dd + 1;
    numTargetDataSACCPrimary(:,dd) = numMostRecentDataSACCPrimary - dd + 1;
    
    % Spd.
    spdSACC{dd} = dataSACC{numTargetDataSACC(dd)}.processedData.P_device;
    spdPrimary1{dd} = dataSACCPrimary1{numTargetDataSACCPrimary(1,dd)}.processedData.P_device;
    spdPrimary2{dd} = dataSACCPrimary2{numTargetDataSACCPrimary(2,dd)}.processedData.P_device;
    spdPrimary3{dd} = dataSACCPrimary3{numTargetDataSACCPrimary(3,dd)}.processedData.P_device;
    
    % Date of measurement.
    dateSACC{dd} = dataSACC{numTargetDataSACC(dd)}.describe.date;
    dateSACCPrimary1{dd} = dataSACCPrimary1{numTargetDataSACCPrimary(1,dd)}.describe.date;
    dateSACCPrimary2{dd} = dataSACCPrimary2{numTargetDataSACCPrimary(2,dd)}.describe.date;
    dateSACCPrimary3{dd} = dataSACCPrimary3{numTargetDataSACCPrimary(3,dd)}.describe.date;
end

%% Plot it.
%
% Set wavelength range.
S = [380 2 201];
wls = SToWls(S);

figure; hold on;
lineColor = {'r--', 'k-'};
lineWidth = 1;
fontSizeTitle = 15;
fontSizeAxis = 12;

% SACC
subplot(4,1,1); hold on;
for dd = 1:nData
    if dd == 1
        lineColor = 'r--';
    else
        lineColor = 'k-';
    end
    plot(wls,spdSACC{dd},lineColor,'linewidth',lineWidth,'DisplayName',num2str(dd));
end
xlabel('Wavelength (nm)','fontsize',fontSizeAxis);
ylabel('Spectral power','fontsize', fontSizeAxis);
title('SACC','fontsize',fontSizeTitle);

% Add legend with date.
%
% Note that when we get a info using get(gca), it sorts the data in a
% reverse order. So, if there is more than one dataset to plot, take care
% to add legend to it.
f = get(gca, 'Children');
idxUpdateInterval = 3;
idxLegend = [1:idxUpdateInterval:idxUpdateInterval*(nData-1)+1];
for ii = 1:nData
    idxContent{ii} = dateSACC{ii};
end
idxContent = flip(idxContent);
legend(f(idxLegend), idxContent);

% SACC primaries
for pp = 1:nPrimaries
    subplot(4,1,pp+1); hold on;
    
    for dd = 1:nData
        if dd == 1
            lineColor = 'r--';
        else
            lineColor = 'k-';
        end
        
        spdTemp = eval(append('spdPrimary',num2str(pp)));
        plot(wls,spdTemp{dd},lineColor,'linewidth',lineWidth,'DisplayName',num2str(dd));
    end
    
    xlabel('Wavelength (nm)','fontsize',fontSizeAxis);
    ylabel('Spectral power','fontsize',fontSizeAxis);
    title(append('SACC Primary ',num2str(pp)),'fontsize',fontSizeTitle);
    
    % Add legend with date.
    f = get(gca, 'Children');
    tempDateSACCPrimary = eval(append('dateSACCPrimary',num2str(pp)));
    
    idxUpdateInterval = 16;
    idxLegend = [1:idxUpdateInterval:idxUpdateInterval*(nData-1)+1];
    for ii = 1:nData
        idxContent{ii} = tempDateSACCPrimary{ii};
    end
    idxContent = flip(idxContent);
    legend(f(idxLegend), idxContent);
end

%% Back to current dir.
cd(dir);
