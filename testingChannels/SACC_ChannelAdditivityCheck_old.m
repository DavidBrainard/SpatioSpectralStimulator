% SACC_ChannelAdditivityCheck
%
% This compares the SPDs between randomly generated spectrum and the sum of
% its single spectrum to check a projecotr additivity.
%
% 10/27/2021 smo   Clean it and makes it more readable.

%% Initialize.
clear; close all;

%% 



%% Load measurement data.
% This is not an elaborate way to read the data, but here it does in this
% way... This part should be updated later.
% One '.mat' file matches a set of spectra.
% To load different file, you need to change the name in 'dataName' and
% 'dataChar' to match the name of the files.
dataName = 'Dataset9)';
dataChar = '(AcrossPrimary)_black10_normal';
load(append(dataName,'Spectrum_single',dataChar,'.mat')); % fw_single
load(append(dataName,'Spectrum_black',dataChar,'.mat')); % fw_blk
load(append(dataName,'Spectrum_white',dataChar,'.mat')); % fw_white
load(append(dataName,'SpectrumInput',dataChar,'.mat')); % SpectrumInput
load(append(dataName,'Spectrum_rand',dataChar,'.mat')) % fw_rand

% Set dataType either within or across primary measurement data.
% dataType = true (within primary) / false (across primary)
dataType = false;

% Set parameters.
S = [380 2 201];
wls = SToWls(S);
nInputLevels = 253;
nPrimaries = 3;
nTest = size(fw_rand,2); % Num of test spectra
subPrimaryWorkingRange = [1:8,10:16];
nChannels = size(subPrimaryWorkingRange,2);

% Delete the not working channel channel.

% Set channel working range for single spectrum and spectrum input for
% random spectra.
if (dataType)
    % Single spectrum
    fw_single = fw_single(:,subPrimaryWorkingRange);
    % Spectrum input
    spectrumInput = SpectrumInput(:,subPrimaryWorkingRange);
else
    % Single spectrum
    subPrimaryWorkingRangeAcrossPrimary = [subPrimaryWorkingRange, subPrimaryWorkingRange+(nChannels+1), subPrimaryWorkingRange+((nChannels+1)*2)];
    fw_single = fw_single(:,subPrimaryWorkingRangeAcrossPrimary);
    fw_singles{1} = fw_single(:,1:nChannels); % Primary 1
    fw_singles{2} = fw_single(:,1+nChannels:2*nChannels); % Primary 2
    fw_singles{3} = fw_single(:,1+2*nChannels:3*nChannels); % Primary 3
    % Spectrum Input
    for pp = 1:nPrimaries
        spectrumInput{pp} = SpectrumInput{pp}(:,subPrimaryWorkingRange);
    end
end

%% Black correction.
% Single spectrum.
if (dataType)
    fw_singleTemp = fw_single - fw_blk;
    fw_singles = max(fw_singleTemp,0);
else
    for pp = 1:nPrimaries
        fw_singlesTemp{pp} = fw_singles{pp} - fw_blk;
        fw_singles{pp} = max(fw_singlesTemp{pp},0);
    end
end

% Random spectrum.
for tt = 1:nTest
    fw_randTemp = fw_rand(:,tt)-fw_blk;
    fw_randTemp = max(fw_randTemp,0);
    fw_rands(:,tt) = fw_randTemp;
end

%% Calculate the sum of single spectra.
fw_blank = zeros(S(3),1);
fw_blanks = zeros(S(3),nChannels);

if (dataType) % within primary
    % Find index which channel is turn-on per each test spectrum.
    idxSpectrumInput = logical(spectrumInput == (nInputLevels-1));
    % Set the spectrum here.
    for tt = 1:nTest
        for ss = 1:nChannels
            if idxSpectrumInput(tt,ss) == 0
                fw_randTemp(:,ss) = fw_blank;
            elseif idxSpectrumInput(tt,ss) == 1
                fw_randTemp(:,ss) = fw_singles(:,ss);
            end
        end
        fw_randSumTemp = sum(fw_randTemp,2);
        fw_randSum(:,tt) = fw_randSumTemp;
    end
    
else (dataType) % across primary
    % Find index which channel is turn-on per each test spectrum.
    for pp = 1:nPrimaries
        idxSpectrumInput{pp} = logical(spectrumInput{pp} == (nInputLevels-1));
    end
 
    fw_randTemp = cell(1,nPrimaries);
    % Set the spectrum here.
    for tt = 1:nTest
        for pp = 1:nPrimaries
            for ss = 1:nChannels
                if idxSpectrumInput{pp}(tt,ss) == 0
                    fw_randTemp{pp}(:,ss) = fw_blank;
                elseif idxSpectrumInput{pp}(tt,ss) == 1
                    fw_randTemp{pp}(:,ss) = fw_singles{pp}(:,ss);
                end
            end
        end
        fw_randSumTemp = sum(cell2mat(fw_randTemp),2);
        fw_randSum(:,tt) = fw_randSumTemp;
    end
end

%% Plot the spectrum.
figure;
for tt = 1:nTest
    subplot(4,nTest/4,tt); hold on;
    plot(wls,fw_rands(:,tt),'k','LineWidth',3)
    plot(wls,fw_randSum(:,tt),'r-','LineWidth',2)
    title(append('Test',num2str(tt)));
    ylim([0 max(max([fw_rand fw_randSum]))]);
end
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
legend('Measurement','Sum result','location','northeast')
title('Comparison of the SPDs between measured and sum result')

%% Compare xy Chromaticity and luminance levels.
% Load color matching function and match the spectrum range.
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,S);

% XYZ calculations.
XYZ_rands = T_xyz * fw_rands;
XYZ_randSum = T_xyz * fw_randSum;
xyY_rands = XYZToxyY(XYZ_rands);
xyY_randSum = XYZToxyY(XYZ_randSum);

% Set color gamut.
colorgamut = XYZToxyY(T_xyz);
colorgamut(:,end+1)=colorgamut(:,1);

% Plot it.
% Luminance
figure; hold on;
testPoints = linspace(1,nTest,nTest);
plot(testPoints,xyY_rands(3,:),'k.','markersize',13);
plot(testPoints,xyY_randSum(3,:),'r+','markersize',12,'linewidth',1);
xlabel('Test point')
ylabel('Luminance (cd/m2)')
title('Luminance level comparison between measurement and sum result')
legend('Measurement','Sum result','location','southeast');

% CIE (x,y) chromaticity
figure; hold on;
plot(xyY_rands(1,:),xyY_rands(2,:),'k.','markersize',13);
plot(xyY_randSum(1,:),xyY_randSum(2,:),'r+','markersize',12,'linewidth',1);
plot(colorgamut(1,:),colorgamut(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
title('Comparison between measured and sum result on the x,y chromaticity')
legend('Measurement','Sum result');
