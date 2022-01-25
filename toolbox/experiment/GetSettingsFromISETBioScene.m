function [primaryFromISETBioGaborCal,settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject,options)
% Get settings in cal format from ISETBio scence.
%
% Syntax:
%    [primaryFromISETBioGaborCal,settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject)
%
% Description:
%    This gets the gabor image settings from the ISETBio scene. As we want
%    to use ISETBio in SACC project, this is for getting back the gabor
%    contrast in a cal format, which should be matched to the contrast
%    gabor before putting inside the ISETBio world.
%
% Inputs: 
%    screenCalObjFromISETBio         - Screen calibration object acquired
%                                      from the ISETBio.
%    ISETBioGaborCalObject           - Structure with the gabor image in a
%                                      cal format that acquired from the
%                                      ISETBio scene.
%    standardGaborCalObject          - Structure with the gabor contrasts
%                              	       and settings in a cal format.
% 
% Outputs:
%    primaryFromISETBioGaborCal      - 
%    settingsFromISETBioGaborCal     - 
%
% Optional key/value pairs:
%    verbose                         - Boolean. Default true. Controls
%                                    plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio, MakeISETBioSceneFromImage

% History:
%   01/21/22  dhb,gka,smo            - Wrote it.
%   01/24/22  smo                    - Made it work.

%% Set parameters.
arguments
    screenCalObjFromISETBio
    ISETBioGaborCalObject
    standardGaborCalObject
    options.verbose (1,1) = true
end

%% Calculate the settings from the ISETBio.
primaryFromISETBioGaborCal = screenCalObjFromISETBio.get('P_device') \ ...
    (ISETBioGaborCalObject.ISETBioGaborCal - screenCalObjFromISETBio.get('P_ambient'));

settingsFromISETBioGaborCal = PrimaryToSettings(screenCalObjFromISETBio,primaryFromISETBioGaborCal);

% Check if it is reasonable.
if (max(abs(standardGaborCalObject.standardSettingsGaborCal(:)-settingsFromISETBioGaborCal(:))./standardGaborCalObject.standardSettingsGaborCal(:)) > 1e-6)
    error('Cannot get home again in settings land');
end

end