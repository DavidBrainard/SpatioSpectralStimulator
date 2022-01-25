function [primaryFromISETBioGaborCal,settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject,options)
% Get settings in cal format from ISETBio scence.
%
% Syntax:
%    [primaryFromISETBioGaborCal,settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject)
%
% Description:
%    TBD
%
% Inputs: 
%    screenCalObjFromISETBio         -
%    ISETBioGaborCalObject           -
%    standardGaborCalObject          -
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
%    SpectralCalISETBio

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

%%
primaryFromISETBioGaborCal = screenCalObjFromISETBio.get('P_device') \ ...
    (ISETBioGaborCalObject.ISETBioGaborCal - screenCalObjFromISETBio.get('P_ambient'));

settingsFromISETBioGaborCal = PrimaryToSettings(screenCalObjFromISETBio,primaryFromISETBioGaborCal);

% Check if it is reasonable.
if (max(abs(standardGaborCalObject.standardSettingsGaborCal(:)-settingsFromISETBioGaborCal(:))./standardGaborCalObject.standardSettingsGaborCal(:)) > 1e-6)
    error('Cannot get home again in settings land');
end

end