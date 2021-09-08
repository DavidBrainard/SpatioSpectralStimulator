function SpatioSpectralStimulatorLocalHook
% ForcedChoiceCMLocalHook
%
% Configure things for working on the SpatioSpectralStimulator project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUseProject('ColorMaterial') to set up for
% this project.  You then edit your local copy to match your configuration.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Define project
projectName = 'SpatioSpectralStimulator';

%% Say hello
fprintf('Running %s local hook\n',projectName);

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify project location
projectBaseDir = tbLocateProject(projectName);

% If we ever needed some user/machine specific preferences, this is one way
% we could do that.
sysInfo = GetComputerInfo();
switch (sysInfo.localHostName)
    case 'eagleray'
        % DHB's desktop
        baseDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
 
    otherwise
        % Some unspecified machine, try user specific customization
        switch(sysInfo.userShortName)
            % Could put user specific things in, but at the moment generic
            % is good enough.
            otherwise % This is for Linux box ('home' instead of 'users') / Semin
                baseDir = fullfile('/home/',sysInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
        end
end

%% Project prefs
setpref(projectName,'LEDSpectraDir',fullfile(baseDir,'SACC_materials','JandJProjector','LEDSpectrumMeasurements'));

% Calibration
setpref('BrainardLabToolbox','CalDataFolder',fullfile(baseDir,'SACC_materials','Calibration'));

%% Set preferences for project output
%
% This will need to be locally configured.
% setpref(projectName,'simulatedDataDir',fullfile(projectBaseDir,'SimulatedData'));
% setpref(projectName,'rayleighDataDir',fullfile(baseDir,'MELA_data','Experiments',projectName,'OLRayleighMatch'));
% setpref(projectName,'rayleighAnalysisDir',fullfile(baseDir,'MELA_analysis','Experiments',projectName, 'OLRayleighMatch'));
% setpref(projectName,'mainExpDir',projectBaseDir);
% setpref(projectName,'analysisDir',fullfile(baseDir,'CNST_analysis',projectName));
% setpref(projectName,'stimulusFolder',fullfile(baseDir,'CNST_materials',projectName,'E3'));
% setpref(projectName,'dataFolder',fullfile(baseDir,'CNST_data',projectName));
% setpref(projectName,'demoDataDir',fullfile(baseDir,'CNST_analysis',projectName,'DemoData'));
% setpref(projectName,'mainCodeDir',fullfile('/Users/', sysInfo.userShortName, 'Documents/MATLAB/projects/Experiments/ColorMaterial/code'));
% setpref(projectName,'calFileName','ColorMaterialCalibration');
% setpref('OneLightToolbox', 'OneLightCalData',fullfile(baseDir,'MELA_materials','Experiments',projectName,'OneLightCalData'));
% setpref('BrainardLabToolbox','CalDataFolder',fullfile(baseDir,'MELA_materials','Experiments',projectName,'OneLightCalData'));
