function [screenCalObj,channelCalObjs] = LoadAndSetExperimentCalFiles(colorDirectionParams,options)
% Load calibration files and refit the gamma.
%
% Syntax:
%    [screenCalObj,channelCalObjs,screenGammaMethod] = LoadAndSetExperimentCalFiles(colorDirectionParams)
%
% Description:
%    It loads the calibration data, both screen and channel, and it does
%    quantized conversion by refitting the gamma.
%
% Inputs:
%    colorDirectionParams         - 
%
% Outputs:
%    screenCalObj                 - 
%    channelCalObjs               - 
%    screenGammaMethod            - 
%
% Optional key/value pairs:
%    options.screenGammaMethod    - 
%    options.channelGammaMethod   - 
%    options.verbose              - 
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,gka,smo         - Wrote it.
%   01/24/22  smo                 - Made it work.

%% Set parameters.
arguments
    colorDirectionParams
    options.screenGammaMethod (1,1) = 2
    options.channelGammaMethod (1,1) = 2
    options.verbose (1,1) = true
end

%% Load calibration data here.
%
% Load screen calibration data.
screenCalObj = LoadCalibration(colorDirectionParams.screenCalName,...
    colorDirectionParams.screenNInputLevels,'setGammaFitMethod',true);

% Load channel calibration data.
nScreenPrimaries = size(colorDirectionParams.channelCalNames,2);
for pp = 1:nScreenPrimaries
    channelCalObjs{pp} = LoadCalibration(colorDirectionParams.channelCalNames{pp},...
        colorDirectionParams.channelNInputLevels,'setGammaFitMethod',false);
end

if (options.verbose)
    disp('Calibration data has been loaded sucessfully!');
end

%% Refit the channel gamma.
%
% Before refitting the gamma, check if the wavelength range matches between
% colorDirectionParams and the calbibration data.
for pp = 1:nScreenPrimaries
    Scheck(pp,:) = channelCalObjs{pp}.get('S');
end
if (any(colorDirectionParams.S ~= Scheck))
    error('Mismatch between calibration file S and that specified at top');
end

% Use quantized conversion from here on.
for cc = 1:3
    SetGammaMethod(channelCalObjs{cc},options.channelGammaMethod);
end

%% Set screen gamma method.
%
% If we set to 0, there is no quantization and the result is excellent. If
% we set to 2, this is quantized at 256 levels and the result is more of a
% mess. The choice of 2 represents what we think will actually happen
% since the real device is quantized.
%
% The point cloud method below reduces this problem.
SetGammaMethod(screenCalObj,options.screenGammaMethod);

end
