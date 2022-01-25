function [calObj] = LoadCalibration(calName,nInputLevels,options)
% Load calibration files and refit its gamma.
%
% Syntax:
%    [calObj] = LoadCalibration(calName,nInputLevels)
%
% Description:
%    Load the calibration measurement files using the designated file names
%    and refit the gamma. You can choose gamma fit method if needed. It can
%    be used for both screen and channel calibrations.
%
% Inputs:
%    calName                 - Calibration file name in string. For SACC
%                              project, it could be either SACC or
%                              SACCPrimary1, SACCPrimary2, SACCPrimary3
%    nInputLevels            - The number of input levels for the
%                              calibration target.
% 
% Outputs:
%    calObj                  - Calibration object in the calibration struct.
%
% Optional key/value pairs:
%    setGammaFitMethod       - Set the gamma fitting method if needed. We
%                              set this for the screen calibration in SACC
%                              project.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/19/22  smo            - Wrote it

%% Set parameters.
arguments
    calName {mustBeMember(calName,{'SACC','SACCPrimary1','SACCPrimary2','SACCPrimary3'})}
    nInputLevels (1,1)
    options.gammaFitMethod = 'identity'
end

%% Load calibration and refit its gamma.
%
% Load the calibration file and make it in cal struct format.
cal = LoadCalFile(calName);
calObj = ObjectToHandleCalOrCalStruct(cal);

% Set the gamma fit method if needed.
if (~isempty(options.gammaFitMethod))
    calObj.set('gamma.fitType',options.gammaFitMethod);
end

% Fit the gamma here.
CalibrateFitGamma(calObj, nInputLevels);

end