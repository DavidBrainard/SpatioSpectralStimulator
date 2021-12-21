% SACC_ChannelAdditivityCheck
%
% This compares the SPDs between randomly generated spectrum and the sum of
% its single spectrum to check a projector additivity both within and
% across screen primary.

% History:
%    10/27/21 smo   Clean it and makes it more readable.
%    12/15/21 smo   Add measurement part. It used to contain only read
%                   the data and plot it.
%    12/20/21 smo   Make it clearer and simpler.

%% Initialize.
clear; close all;

%% Set parameters here.
%
% This part would stay the same.
S = [380 2 201];
wls = SToWls(S);
nInputLevels = 253;
nPrimaries = 3;
nChannels = 16;

% This part can be changable.
channelIntensity = 1;
arbitraryBlack = 0.05;
nTestSamples = 20;
nTestSamplePeaks = 3;

% Set the test type ('within' or 'across' screen primary).
TESTTYPE = 'within';

%% Open the screen and connect spectroradiometer.
%
% The DLP screen will be controlled as white plain screen.
OpenPlainScreen([1 1 1]);
OpenSpectroradiometer;

% Measure an dark ambient setting for black correction.
channelSettingsBlack = ones(nChannels,nPrimaries) * arbitraryBlack;
SetChannelSettings(channelSettingsBlack);
darkAmbient = MeasureSpectroradiometer;

switch TESTTYPE
    case 'within'
        %% Measure single spectrum.
        for cc = 1:nChannels
            channelSettingsSingle = ones(nChannels,nPrimaries) .* arbitraryBlack;
            channelSettingsSingle(cc,:) = channelIntensity;
            SetChannelSettings(channelSettingsSingle);
            
            % Measure it.
            spdSingle(:,cc) = MeasureSpectroradiometer;
            fprintf('Single channel measurement complete! - (%d/%d)\n',cc,nChannels);
            
            % Black correction here.
            spdSingle(:,cc) = spdSingle(:,cc) - darkAmbient;
        end
        
        %% Generate and measure test random spectrum.
        %
        % Make index for creating random spectrum.
        for ii = 1:nTestSamples
            % Randomize the channels to turn on.
            idxTestSamplePeaks(:,ii) = randsample([1:1:nChannels],nTestSamplePeaks);
        end
        
        % Measure it.
        for ii = 1:nTestSamples
            for pp = 1:nTestSamplePeaks
                channelSettingsTestSample = ones(nChannels, nPrimaries) * arbitraryBlack;
                channelSettingsTestSample(idxTestSamplePeaks(pp,ii)) = channelIntensity;
            end
            SetChannelSettings(channelSettingsTestSample);
            spdTestSample(:,ii) = MeasureSpectroradiometer;
            fprintf('Test sample measurement complete! - (%d/%d)\n',ii,nTestSamples);
        end
        
        % Calculate the sum of spectra matching the combinations of test samples.
        for ii = 1:nTestSamples
            spdTestSampleSum(:,ii) = sum(spdSingle(idxTestSamplePeaks(:,ii)));
        end
        
        %% XYZ calculations.
        %
        % Load color matching function and match the spectrum range.
        load T_xyzJuddVos;
        T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,S);
        
        % XYZ calculations.
        XYZTestSample = T_xyz * spdTestSample;
        XYZTestSampleSum = T_xyz * spdTestSampleSum;
        xyYTestSample = XYZToxyY(XYZTestSample);
        xyYTestSampleSum = XYZToxyY(XYZTestSampleSum);
        
        % Set spectral locus.
        spectralLocus = XYZToxyY(T_xyz);
        spectralLocus(:,end+1) = spectralLocus(:,1);
        
    case 'across'
        %% Measure single spectrum.
        for pp = 1:nPrimaries
            for cc = 1:nChannels
                % Set target channel setting and the others as arbitrary black.
                % We will correct the black later on.
                channelSettingsSingle = ones(nChannels,nPrimaries) * arbitraryBlack;
                channelSettingsSingle(cc,pp) = channelIntensity;
                SetChannelSettings(channelSettingsSingle);
                
                % Measure it.
                spdSingle(pp,cc,:) = MeasureSpectroradiometer;
                fprintf('Single channel measurement complete! - screen primary:%d / channel:(%d/%d)\n',pp,cc,nChannels);
                
                % Black correction here.
                spdSingle(pp,cc,:) = spdSingle(pp,cc,:) - darkAmbient;
            end
        end
        
        %% Generate and measure test random spectrum.
        %
        % Make index for creating random spectrum.
        for ii = 1:nTestSamples
            % Randomize the channels to turn on.
            idxTestSamplePeaks(:,ii) = randsample([1:1:nChannels],nTestSamplePeaks);
        end
        
        % Measure it.
        for ii = 1:nTestSamples
            channelSettingsTestSample = ones(nChannels, nPrimaries) * arbitraryBlack;
            
            % Set the peak here.
            for ss = 1:nPrimaries
                channelSettingsTestSample(idxTestSamplePeaks(pp,ii),ss) = channelIntensity;
            end
            
            SetChannelSettings(channelSettingsTestSample);
            spdTestSample(:,ii) = MeasureSpectroradiometer;
            fprintf('Test sample measurement complete! - (%d/%d)\n',ii,nTestSamples);
        end
        
        % Calculate the sum of spectra matching the combinations of test samples.
        for ii = 1:nTestSamples
            for ss = 1:nPrimaries
                spdTestSampleSumTemp = spdSingle(ss, idxTestSamplePeaks(ss,ii), :);
                spdTestSampleSum(:,ii) = sum(spdTestSampleSumTemp);
            end
        end
    otherwise
end

%% Plot it.
figure;
for tt = 1:nTestSamples
    subplot(5,nTestSamples/5,tt); hold on;
    plot(wls,spdTestSample(:,tt),'k','LineWidth',3)
    plot(wls,spdTestSampleSum(:,tt),'r-','LineWidth',2)
    title(append('Test',num2str(tt)));
    ylim([0 max(max([spdTestSample spdTestSampleSum]))]);
end
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
legend('Measurement','Sum result','location','northeast')
title('Comparison of the SPDs between measured and sum result')

% Luminance.
figure; hold on;
numTestSamples = linspace(1,nTest,nTest);
plot(numTestSamples,xyYTestSample(3,:),'k.','markersize',13);
plot(numTestSamples,xyYTestSampleSum(3,:),'r+','markersize',12,'linewidth',1);
xlabel('Test point')
ylabel('Luminance (cd/m2)')
title('Luminance level comparison between measurement and sum result')
legend('Measurement','Sum result','location','southeast');

% CIE (x,y) chromaticity
figure; hold on;
plot(xyYTestSample(1,:),xyYTestSample(2,:),'k.','markersize',13);
plot(xyYTestSampleSum(1,:),xyYTestSampleSum(2,:),'r+','markersize',12,'linewidth',1);
plot(spectralLocus(1,:),spectralLocus(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
title('Comparison between measured and sum result on the x,y chromaticity')
legend('Measurement','Sum result');

%% Save the data.
%
% Close the screen and projector.
CloseScreen;
CloseSpectroradiometer;

% Save data with the name containing dayTimestr.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    testFilename = fullfile(testFiledir,sprintf('additivityCheck_%s_%s',TESTTYPE,dayTimestr));
    save(testFilename,'spdSingle','spdTestSample','spdTestSampleSum','nTestSamples','nTestSamplePeaks','idxTestSamplePeaks','S');
end
