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

%% Initialize.
clear; close all;

%% Set variables.
frequencyFlicker = 25;
leftButton = true;
gaussianWindow = false;
bgColor = 'white';
SAVETHERESULTS = true;

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
        data = GetMatchingRedForRGFlicker('calibrate',true,'bgColor',bgColor,...
            'gaussianWindow',gaussianWindow);
        
    case 'demo'
        % 1) Look at stimulus demo.
        %
        % For practice, subject will see the flicker stimulus by controlling the
        % red starting either from top or bottom.
        nTrials = 2;
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
        dataMain = GetMatchingRedForRGFlicker('nTrials',nTrials,'bgColor',bgColor,...
            'leftButton',leftButton,'gaussianWindow',gaussianWindow,'frequencyFlicker',frequencyFlicker);
end

%% Collect all data and save it.
if (SAVETHERESULTS)
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        
        % Set the file name and save.
        dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
        testFilename = fullfile(testFiledir,subjectName,...
            sprintf('Flicker_%s_%s',subjectName,dayTimestr));
        save(testFilename,'dataMain');
    end
end
