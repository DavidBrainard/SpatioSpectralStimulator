function [] = GetSettingsFromISETBioScene()
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
end

%%

primaryFromISETBioGaborCal = screenCalObjFromISETBio.get('P_device') \ ...
    (ISETBioGaborCal-screenCalObjFromISETBio.get('P_ambient'));
settingsFromISETBioGaborCal = PrimaryToSettings(screenCalObjFromISETBio,primaryFromISETBioGaborCal);
if (max(abs(standardSettingsGaborCal(:)-settingsFromISETBioGaborCal(:))./standardSettingsGaborCal(:)) > 1e-6)
    error('Cannot get home again in settings land');
end

end