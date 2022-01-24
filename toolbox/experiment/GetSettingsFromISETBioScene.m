function [primaryFromISETBioGaborCal,settingsFromISETBioGaborCal] = GetSettingsFromISETBioScene(screenCalObjFromISETBio,ISETBioGaborCalObject,standardGaborCalObject)
% d
%
% Syntax:
%    d
%
% Description:
%    d
%
% Inputs:
%    d                       -
%
% Outputs:
%    d                       -
%
% Optional key/value pairs:
%    d                       - d
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,ga,smo     - Wrote it

%% Set parameters.
arguments
    screenCalObjFromISETBio
    ISETBioGaborCalObject
    standardGaborCalObject
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