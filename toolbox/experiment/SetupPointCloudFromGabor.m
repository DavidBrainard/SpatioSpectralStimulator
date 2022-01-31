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
%   01/31/22  smo                   - It is possible to make multiple
%                                     target contrast gabors inside this
%                                     function.

%% Set parameters.
arguments
    colorDirectionParams
    rawMonochromeContrastGaborCal
    screenCalObj
    screenBgExcitations
    options.verbose (1,1) = true
end

%% Set up point cloud of contrasts for all possible settings.
%
% Make a point cloud here. It will take a while.
if (options.verbose)
    disp('Starting to make contrast point cloud...');
end
[contrastPtCld, ptCldSettingsCal] = SetupContrastPointCloud(screenCalObj,screenBgExcitations,'verbose',options.verbose);
% Save in a struct.
ptCldObject.contrastPtCld = contrastPtCld;
ptCldObject.ptCldSettingsCal = ptCldSettingsCal;
if (options.verbose)
    disp('Contrast point cloud has been successfully created!');
end

%% Get cone contrast/excitation gabor image.
%
% Scale target cone contrast vector at max excursion by contrast modulation
% at each pixel.  This is done by a single matrix multiply plus a lead
% factor.  We work cal format here as that makes color transforms
% efficient.
nContrastPoints = size(colorDirectionParams.spatialGaborTargetContrast,2);

% Make a loop for the number of target gabor contrasts.
for cc = 1:nContrastPoints
    desiredContrastGaborCal = colorDirectionParams.spatialGaborTargetContrast(cc) * colorDirectionParams.targetStimulusContrastDir * rawMonochromeContrastGaborCal;
    
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
    
    % Save the results in a struct.
    standardGaborCalObject.desiredContrastGaborCal{cc} = desiredContrastGaborCal;
    standardGaborCalObject.desiredExcitationsGaborCal{cc} = desiredExcitationsGaborCal;
    standardGaborCalObject.standardPrimariesGaborCal{cc} = standardPrimariesGaborCal;
    standardGaborCalObject.desiredSpdGaborCal{cc} = desiredSpdGaborCal;
    standardGaborCalObject.standardSettingsGaborCal{cc} = standardSettingsGaborCal;
    standardGaborCalObject.standardPredictedPrimariesGaborCal{cc} = standardPredictedPrimariesGaborCal;
    standardGaborCalObject.standardPredictedExcitationsGaborCal{cc} = standardPredictedExcitationsGaborCal;
    standardGaborCalObject.standardPredictedContrastGaborCal{cc} = standardPredictedContrastGaborCal;
end

end
