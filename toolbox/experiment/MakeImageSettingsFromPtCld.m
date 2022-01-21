function [] = MakeImageSettingsFromPtCld()
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

%% We want this routine to take contrast explicitly, expressed relative to max contrast we set up,
% when it makes the image.  We will call this multiple times to make
% stimuli of different contrasts.

%% Get image from point cloud in cal format.
uniqueQuantizedSettingsGaborCal = SettingsFromPointCloud(contrastPtCld,desiredContrastGaborCal,ptCldSettingsCal);

% Print out min/max of settings
fprintf('Gabor image min/max settings: %0.3f, %0.3f\n',min(uniqueQuantizedSettingsGaborCal(:)), max(uniqueQuantizedSettingsGaborCal(:)));

% Get contrasts we think we have obtianed
uniqueQuantizedExcitationsGaborCal = SettingsToSensor(screenCalObj,uniqueQuantizedSettingsGaborCal);
uniqueQuantizedContrastGaborCal = ExcitationsToContrast(uniqueQuantizedExcitationsGaborCal,screenBgExcitations);

% Plot of how well point cloud method does in obtaining desired contrats
figure; clf;
plot(desiredContrastGaborCal(:),uniqueQuantizedContrastGaborCal(:),'r+');
axis('square');
xlabel('Desired L, M or S contrast');
ylabel('Predicted L, M, or S contrast');
title('Quantized unique point cloud image method');

%% Convert representations we want to take forward to image format
desiredContrastGaborImage = CalFormatToImage(desiredContrastGaborCal,stimulusN,stimulusN);
standardPredictedContrastImage = CalFormatToImage(standardPredictedContrastGaborCal,stimulusN,stimulusN);
standardSettingsGaborImage = CalFormatToImage(standardSettingsGaborCal,stimulusN,stimulusN);
uniqueQuantizedContrastGaborImage = CalFormatToImage(uniqueQuantizedContrastGaborCal,stimulusN,stimulusN);
end