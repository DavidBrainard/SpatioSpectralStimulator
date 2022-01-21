function [] = LoadAndSetExperimentCalFiles()
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

%% Load screen calibration and refit its gamma.
%
% Load screen calibration.
screenCalObj = LoadCalibration(colorDirectionParams.screenCalName,...
    colorDirectionParams.screenNInputLevels,'setGammaFitMethod',true);

% Load channel calibration.
nScreenPrimaries = size(colorDirectionParams.channelCalNames,2);
for pp = 1:nScreenPrimaries
    channelCalObjs{pp} = LoadCalibration(colorDirectionParams.channelCalNames{pp},...
        colorDirectionParams.channelNInputLevels,'setGammaFitMethod',false);
end

%% Get out some data to work with.
%
% This is from the channel calibration file.
for pp = 1:nScreenPrimaries
    Scheck(pp,:) = channelCalObjs{pp}.get('S');
end
if (any(colorDirectionParams.S ~= Scheck))
    error('Mismatch between calibration file S and that specified at top');
end
wls = SToWls(colorDirectionParams.S);
nChannels = channelCalObjs{1}.get('nDevices');

%% Use quantized conversion from here on.
%
% Comment in the line that refits the gamma to see
% effects of extreme quantization one what follows.
%
% CalibrateFitGamma(channelCalObjs{1},10);
channelGammaMethod = 2;
for cc = 1:3
    SetGammaMethod(channelCalObjs{cc},channelGammaMethod);
end

%% Set screen gamma method.
%
% If we set to 0, there is no quantization and the result is excellent.
% If we set to 2, this is quantized at 256 levels and the result is more
% of a mess.  The choice of 2 represents what we think will actually happen
% since the real device is quantized.
%
% The point cloud method below reduces this problem.
screenGammaMethod = 2;
SetGammaMethod(screenCalObj,screenGammaMethod);
end