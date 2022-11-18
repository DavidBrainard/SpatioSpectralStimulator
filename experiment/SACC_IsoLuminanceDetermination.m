% SACC_IsoLuminanceDetermination
%
% This is for iso-luminance determination for SACC project. Here we display
% simple flicker (red-green) stimuli for observers to adjust the intensity
% of red until the flicker disappears.

% History:
%    08/11/22   smo        - Started on it.
%    08/17/22   smo        - Added a control part to update the intensity of
%                            red light while making red-green flicker.
%    08/18/22   dhb, smo   - Now code is working fine.
%    08/19/22   smo        - Added Gaussian noise on the white background.
%    09/14/22   dhb, smo   - Flicker frequency (frame numbers) has been
%                            corrected.
%    11/08/22   smo        - Flicker session has been substituted with the
%                            function.
%    11/18/22   smo        - Added a calibartion mode.

%% Initialize.
clear; close all;

%% Set variables.
frequencyFlicker = 25;
leftButton = true;
gaussianWindow = false;
bgColor = 'white';
SAVETHERESULTS = true;
verbose = true;

%% Get subject name and which mode to run.
%
% Which mode to run.
while 1
    inputMessageMode = 'Which mode to run [calibrate, demo, practice, main]: ';
    whichMode = input(inputMessageMode, 's');
    whichModeOptions = {'calibrate', 'demo', 'practice', 'main'};
    
    if ismember(whichMode, whichModeOptions)
        fprintf('\t Flicker code will be run in (%s) mode! \n', whichMode);
        break
    end
    
    disp('Running mode should be selected within [demo, practice, main]!');
end

% Subject name.
if strcmp(whichMode, 'main')
    inputMessageName = 'Enter subject name: ';
    subjectName = input(inputMessageName, 's');
elseif strcmp(whichMode, 'calibrate')
    subjectName = 'calibrate';
end

%% Run the flicker code according to different mode.
switch whichMode
    case 'calibrate'
        calibrate = true;
        
        % Measurement happens here.
        data = GetMatchingRedForRGFlicker('calibrate',true,'bgColor',bgColor,...
            'gaussianWindow',gaussianWindow);
        
        % Black correction.
        spdGreenBC = data.spdGreen - data.spdAmbient;
        spdRedBC = data.spdRed - data.spdAmbient;
        
        % Set wavelength range and calculate XYZ.
        S = [380 2 201];
        wls = SToWls(S);

        load T_xyzJuddVos;
        T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,S);

        xyzRedBC = T_xyz*spdRedBC;
        xyzGreenBC = T_xyz*spdGreenBC;
        
        % Plot it.
        if (verbose)
            % Spectrum.
            figure; clf; hold on;
            plot(wls, spdGreenBC, 'g--', 'linewidth', 2);
            plot(wls, spdRedBC, 'r-', 'linewidth', 0.5);
            xlabel('Wavelength (nm)','fontsize',13)
            ylabel('Spectral power','fontsize',13);
            legend('Green', 'Red','fontsize',13);
            xlim([380 780]);
            xticks([380:80:780]);

            % Settings vs. Output.
            figure; clf; hold on;
            plot(data.greenIntensity, xyzGreenBC(2),'ko','markersize',11,'markerfacecolor','g');
            plot(data.redIntensity, xyzRedBC(2,:),'ko','markersize',11,'markerfacecolor','r')
            xlabel('Input settings (8-bit)','fontsize',13)
            ylabel('Output (no unit)','fontsize',13);
            legend('Green', 'Red','fontsize',13,'location','southeast');
        end

    case 'demo'
        % 1) Look at stimulus demo.
        %
        % For practice, subject will see the flicker stimulus by controlling the
        % red starting either from top or bottom.
        nTrials = 1;
        GetMatchingRedForRGFlicker('nTrials',nTrials,'bgColor',bgColor,...
            'leftButton',leftButton,'gaussianWindow',gaussianWindow,'frequencyFlicker',frequencyFlicker);
        
        % Stop after trials.
        return;
    
    case 'practice'
        % 2) Practice trials.
        %
        % We will repeat the whole session twice so that subject will do a total of
        % 4 flicker sessions (2 repeatitions x 2 starting points either top or bottom).
        nTrials = 2;
        dataPractice = GetMatchingRedForRGFlicker('nTrials',nTrials,'bgColor',bgColor,...
            'leftButton',leftButton,'gaussianWindow',gaussianWindow,'frequencyFlicker',frequencyFlicker);
    
        % Stop after trials.
        return;
    
    case 'main'
        % 3) Main session.
        %
        % By protocol, we will meausre a total of 6 sessions (3 repeatitions x 2
        % starting points either top or bottom).
        nTrials = 6;
        data = GetMatchingRedForRGFlicker('nTrials',nTrials,'bgColor',bgColor,...
            'leftButton',leftButton,'gaussianWindow',gaussianWindow,'frequencyFlicker',frequencyFlicker);
end

%% Collect all data and save it.
if (SAVETHERESULTS)
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        
        % Set the file name and save.
        switch whichMode
            case 'calibrate'
                testFilename = fullfile(testFiledir,'CheckFlickerPhotom',...
                    sprintf('checkFlickerPhotom_%s',dayTimestr));
                save(testFilename,'data');

                SAVECALIBRATIONPLOT = true;
                if (SAVECALIBRATIONPLOT)
                    if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
                        testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),...
                            'CheckFlickerPhotom');
                        testFilename = fullfile(testFiledir,sprintf('CheckFlickerPhotom_%s',dayTimestr));
                        testFileFormat = '.tiff';
                        saveas(gcf,append(testFilename,testFileFormat));
                        fpintf('Figure has been succesfully saved! \n');
                    end
                end
                
            case 'main'
                testFilename = fullfile(testFiledir,subjectName,...
                    sprintf('flickerPhotom_%s_%s',subjectName,dayTimestr));
                save(testFilename,'data');
        end
    end
end
fpintf('Data file (%s) has been succesfully saved! \n', whichMode);
