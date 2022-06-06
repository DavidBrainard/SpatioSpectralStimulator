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
%    01/07/22 smo   Added the option skipping the measurement.

%% Initialize.
clear all; close all;

%% Set parameters here.
%
% This part would stay the same.
S = [380 2 201];
wls = SToWls(S);
nInputLevels = 253;
nPrimaries = 3;
nChannels = 16;

% This part can be changable.
% The number of peaks on each test sample (nTestSamplePeaks) can be changed
% only for within primary test for now. It should be set as 3 when you test
% across primary matching the number of screen primaries.
channelIntensity = 1;
arbitraryBlack = 0.05;
nTestSamples = 20;
nTestSamplePeaks = 3;

% Set the test type ('within' or 'across' screen primary).
TESTTYPE = 'across';
MEASURE = false;

%% Open the screen and connect to spectroradiometer.
%
% The DLP screen will be controlled as white plain screen.
if (MEASURE)
    OpenPlainScreen([1 1 1]);
    OpenSpectroradiometer;

    %% Measure an dark ambient setting for black correction.
    channelSettingsBlack = ones(nChannels,nPrimaries) * arbitraryBlack;
    SetChannelSettings(channelSettingsBlack);
    darkAmbient = MeasureSpectroradiometer;
    disp('Dark ambient measurement complete!');

    %% Further measurement happens from here.
    %
    % Set the TESTTYPE either 'within' or 'across' from the above.
    % Both tests have the same routine which starts from measuring a single
    % spectrum, then generate random spectrum to compare it with the sum of
    % single spectra.
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
                channelSettingsTestSample = ones(nChannels, nPrimaries) * arbitraryBlack;
                for kk = 1:nTestSamplePeaks
                    channelSettingsTestSample(idxTestSamplePeaks(kk,ii),:) = channelIntensity;
                end
                SetChannelSettings(channelSettingsTestSample);
                spdTestSample(:,ii) = MeasureSpectroradiometer;
                fprintf('Test sample measurement complete! - (%d/%d)\n',ii,nTestSamples);

                % Black correction here.
                spdTestSample(:,ii) = spdTestSample(:,ii) - darkAmbient;
            end

            % Calculate the sum of spectra matching the combinations of test samples.
            for ii = 1:nTestSamples
                for pp = 1:nTestSamplePeaks
                    spdTestSampleSumTemp(:,pp) = spdSingle(:, idxTestSamplePeaks(pp,ii));
                end
                spdTestSampleSum(:,ii) = sum(spdTestSampleSumTemp,2);
            end

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
                    spdSingleTemp = MeasureSpectroradiometer;
                    fprintf('Single channel measurement complete! - screen primary:%d / channel:(%d/%d)\n',pp,cc,nChannels);

                    % Black correction here.
                    spdSingle(pp,cc,:) = spdSingleTemp - darkAmbient;
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
                for pp = 1:nPrimaries
                    channelSettingsTestSample(idxTestSamplePeaks(pp,ii),pp) = channelIntensity;
                end

                SetChannelSettings(channelSettingsTestSample);
                spdTestSample(:,ii) = MeasureSpectroradiometer;
                fprintf('Test sample measurement complete! - (%d/%d)\n',ii,nTestSamples);

                % Black correction here.
                spdTestSample(:,ii) = spdTestSample(:,ii) - darkAmbient;
            end

            % Calculate the sum of spectra matching the combinations of test samples.
            for ii = 1:nTestSamples
                for pp = 1:nPrimaries
                    spdTestSampleSumTemp(:,pp) = spdSingle(pp, idxTestSamplePeaks(pp,ii), :);
                end
                spdTestSampleSum(:,ii) = sum(spdTestSampleSumTemp,2);
            end
        otherwise
    end
else
    % Load the data if the measurement is skipped.
    if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
        testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
        testFilename = GetMostRecentFileName(testFiledir,sprintf('additivityCheck_%s',TESTTYPE),'olderDate',0);
        load(testFilename);
    end
end

%% XYZ calculations.
%
% Before calculation of XYZ, round any negative parts in the spectrum
% because of the black correction.
spdSingle(spdSingle < 0) = 0;
spdTestSample(spdTestSample < 0) = 0;
spdTestSampleSum(spdTestSampleSum < 0) = 0;

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

%% Plot it.
%
% All test colors spectra.
figure;
for tt = 1:nTestSamples
    subplot(5,nTestSamples/5,tt); hold on;
    plot(wls,spdTestSample(:,tt),'k-','LineWidth',2)
    plot(wls,spdTestSampleSum(:,tt),'r--','LineWidth',2)
    title(append('Test',num2str(tt)));
    ylim([0 max(max([spdTestSample spdTestSampleSum]))]);
end

% One spectrum.
numTestSpectrum = 1;
figure; hold on;
plot(wls,spdTestSample(:,numTestSpectrum),'k-','LineWidth',2)
plot(wls,spdTestSampleSum(:,numTestSpectrum),'r--','LineWidth',2)
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance')
legend('Measurement','Sum result','location','northeast')

% Luminance.
figure; hold on;
numTestSamples = linspace(1,nTestSamples,nTestSamples);
axisLimit = round(max([xyYTestSample(3,:) xyYTestSampleSum(3,:)], [], 'all'))+1;
plot(xyYTestSample(3,:), xyYTestSampleSum(3,:),'r.','markersize',15);
plot([0 axisLimit],[0 axisLimit],'k-');
xlabel('Luminance - Measurement (cd/m2)')
ylabel('Luminance - Sum result (cd/m2)')
ylim([0 axisLimit]);
xlim([0 axisLimit]);
yticks([0:1:axisLimit]);
xticks([0:1:axisLimit]);
axis('square');
% title('Luminance level comparison between measurement and sum result')
legend('Test spectra','45-deg line','location','southeast');

% CIE (x,y) chromaticity
figure; hold on;
plot(xyYTestSample(1,:),xyYTestSample(2,:),'k.','markersize',13);
plot(xyYTestSampleSum(1,:),xyYTestSampleSum(2,:),'r+','markersize',12,'linewidth',1);
plot(spectralLocus(1,:),spectralLocus(2,:),'k-');
xlabel('CIE x')
ylabel('CIE y')
xlim([0 1]);
ylim([0 1]);
yticks([0:0.2:1]);
xticks([0:0.2:1]);
axis('square');
% title('Comparison between measured and sum result on the x,y chromaticity')
legend('Measurement','Sum result');

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
        testFilename = fullfile(testFiledir,sprintf('additivityCheck_%s_%s',TESTTYPE,dayTimestr));
        save(testFilename,'spdSingle','spdTestSample','spdTestSampleSum','nTestSamples','nTestSamplePeaks',...
            'idxTestSamplePeaks','arbitraryBlack','channelIntensity','S');
    end
end
