function [ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,screenCalObj,screenBgExcitations,options)
% Set up point cloud from the gabor image.
%
% Syntax:
%    [ptCldObject,standardGaborCalObject] = SetupPointCloudFromGabor(colorDirectionParams,rawMonochromeContrastGaborCal,screenCalObj,screenBgExcitations)
%
% Description:
%    This creates the point cloud that has the all possible combinations
%    of the contrasts. This will be used to create a gabor image with a
%    desired contrast.
%
% Inputs:
%    colorDirectionParams           - Structure with the parameters to
%                                     calculate a contrast gabor image.
%    rawMonochromeContrastGaborCal  - Monochrome contrast gabor image in a
%                                     cal format to compute the gabor image
%                                     with desired contrast and color
%                                     direction.
%    screenCalObj                   - Screen calibration object.
%    screenBgExcitations            - Screen background cone excitations.
%
% Outputs:
%    ptCldObject                    - Structure with the contrasts for all
%                                     possible settings using the point
%                                     cloud method.
%    standardGaborCalObject         - Structure with the gabor contrasts
%                                     and settings in a cal format.
%
% Optional key/value pairs:
%    verbose                        - Boolean. Default true. Controls
%                                     plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,gka,smo           - Wrote it.
%   01/24/22  smo                   - Made it work.

%% Set parameters.
arguments
    colorDirectionParams
    rawMonochromeContrastGaborCal
    screenCalObj
    screenBgExcitations
    options.verbose (1,1) = true
end

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

%% Set up point cloud of contrasts for all possible settings.
[contrastPtCld, ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',options.verbose);
if (options.verbose)
    disp('Contrast point cloud has been successfully created!');
end

%% Save the results in a struct to print out.
%
% Note that 'desiredContrastGaborCal' is saved in both objects for
% convenience. It can be changed later on.
%
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
