function [ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,screenCalObj,screenBgExcitations,options)
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
    colorDirectionParams
    rawMonochromeContrastGaborCal
    screenCalObj
    screenBgExcitations
    options.verbose (1,1) = true
end

%%

%% Get cone contrast/excitation gabor image.
%
% Scale target cone contrast vector at max excursion by contrast modulation
% at each pixel.  This is done by a single matrix multiply plus a lead
% factor.  We work cal format here as that makes color transforms
% efficient.
desiredContrastGaborCal = colorDirectionParams.spatialGaborTargetContrast * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastGaborCal;

% Convert cone contrast to excitations
desiredExcitationsGaborCal = ContrastToExcitation(desiredContrastGaborCal,screenBgExcitations);

% Get primaries using standard calibration code, and desired spd without
% quantizing.
standardPrimariesGaborCal = SensorToPrimary(screenCalObj,desiredExcitationsGaborCal);
desiredSpdGaborCal = PrimaryToSpd(screenCalObj,standardPrimariesGaborCal);

% Gamma correct and quantize (if gamma method set to 2 above; with gamma
% method set to zero there is no quantization).  Then convert back from
% the gamma corrected settings.
standardSettingsGaborCal = PrimaryToSettings(screenCalObj,standardPrimariesGaborCal);
standardPredictedPrimariesGaborCal = SettingsToPrimary(screenCalObj,standardSettingsGaborCal);
standardPredictedExcitationsGaborCal = PrimaryToSensor(screenCalObj,standardPredictedPrimariesGaborCal);
standardPredictedContrastGaborCal = ExcitationsToContrast(standardPredictedExcitationsGaborCal,screenBgExcitations);

%% Set up point cloud of contrasts for all possible settings
[contrastPtCld, ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',options.verbose);

%% Save the results in a struct.

% Standard gabor cal object.
standardGaborCalObject.desiredContrastGaborCal = desiredContrastGaborCal;
standardGaborCalObject.desiredExcitationsGaborCal = desiredExcitationsGaborCal;
standardGaborCalObject.standardPrimariesGaborCal = standardPrimariesGaborCal;
standardGaborCalObject.desiredSpdGaborCal = desiredSpdGaborCal;
standardGaborCalObject.standardSettingsGaborCal = standardSettingsGaborCal;
standardGaborCalObject.standardPredictedPrimariesGaborCal = standardPredictedPrimariesGaborCal;
standardGaborCalObject.standardPredictedExcitationsGaborCal = standardPredictedExcitationsGaborCal;
standardGaborCalObject.standardPredictedContrastGaborCal = standardPredictedContrastGaborCal;

% Point cloud object.
ptCldObject.desiredContrastGaborCal = desiredContrastGaborCal;
ptCldObject.contrastPtCld = contrastPtCld;
ptCldObject.ptCldSettingsCal = ptCldSettingsCal;

end
