function [isolatingPrimaries,isolatingPrimariesQuantized,isolatingSpd,isolatingContrast,isolatingModulationPrimaries] = FindChannelPrimaries(targetMaxLMSContrast,targetPrimaryHeadroom,targetContrastReMax,bgPrimaries, ...
    T_cones,subprimaryCalstructData,B_natural,projectIndices,primaryHeadroom,targetLambda,options)
% Find the desired target primaries so that the three primaries can produce
% desired LMS contrast.
%
% Syntax:
%    [isolatingPrimaries,isolatingPrimariesQuantized,isolatingSpd,isolatingContrast] = FindChannelPrimaries(targetMaxLMSContrast,targetPrimaryHeadroom,targetContrastReMax,bgPrimaries, ...
%                                                                                      T_cones,subprimaryCalstructData,B_natural,projectIndices,primaryHeadroom,targetLambda,options)
%
% Description:
%    This function searches for three primaries based on the desired LMS
%    contrast from the background. The shape of the final obtained spectrums will
%    depend on the variable 'targetLambda' which constrains the smoothness
%    of the curve.
%
%    This routine allows for optional passing of an ExtraAmbientSpd
%    variable, to account for use when we are calculating for subprimaries
%    for one primary channel, knowing that the other two will also
%    contribute to the overall background.
%
% Inputs:
%    targetMaxLMSContrast -       The maximum LMS contrast values for the target
%                                 are stored here for all three primaries.
%                                 Each column should be aligned as each
%                                 primary.
%    targetPrimaryHeadroom -      This decides a little headroom when setting
%                                 up the target LMS contrast to avoid a
%                                 numerical error at the edges.
%    targetContrastReMax -        The amount of the contrast that we want
%                                 to use (we do not want the whole range).
%    bgPrimaries -                The primaries of the background.
%    T_cones -                    Spectral cone sensitivities in standard PTB format.
%    subPrimaryCalstructData -    The subprimary calstructdata to calculate
%                                 some routines including 'PrimaryToSpd','SettingsToPrimary',
%                                 and 'PrimaryToSettings'
%    B_natural -                  Basis functions to enforce the smoothness
%                                 on the desired background spectrum.
%                                 This also should be in standard PTB
%                                 format.
%    projectIndices -             Define wavelength range that will be used to enforce the smoothnes
%                                 thorugh the projection onto an underlying basis set.
%    primaryHeadroom -            Set the headroom to prevent running into
%                                 numerical error at the edges.
%    targetLambda -               This number decides the level of smoothness of
%                                 the final calculation result of the
%                                 desired background spectrum. Higher value
%                                 makes the curve more smoother,
%                                 but it can possibly lose accuracy
%                                 in matching the desired xy chromaticity.
%                                 You may want to play with targetLambda a bit
%                                 by hand to determine tradeoff between
%                                 smoothness and obtaining desired chromaticity.
%
% Outputs:
%    isolatingPrimaries -         Obtained isolating primaries which
%                                 reproduce the desired LMS contrasts.
%    isolatingPrimariesQuantized -The quantized results of the isolating primareis.
%    isolatingSpd -               The spectrums of the isolating primaries.
%                                 This is calculated based on the quantized
%                                 primaries results ('isolatingPrimariesQuantized').
%    isolatingContrast -          The LMS contrasts that isolating primaries
%                                 reproduce. This is also calculated
%                                 based on the quantized isolating
%                                 primaries.
%    isolatingModulationPrimaries - The differential primaries that head in
%                                 the desired direction.
%
% Optional key/value pairs:
%    'ExtraAmbientSpd' -          Vector specifying extra ambient light spd
%                                 that needs tp be accounted for, for
%                                 example from the background from other primary channels.
%
% History:
%    09/29/21  smo                Started on it
%    10/05/21  dhb,smo            Clean it and added the feature of 'ExtraAmbientSpd'.
%    10/14/21  dhb                Quantize computation of background primaries before computing spectra.

%% Set parameters.
arguments
    targetMaxLMSContrast
    targetPrimaryHeadroom
    targetContrastReMax
    bgPrimaries
    T_cones
    subprimaryCalstructData
    B_natural
    projectIndices
    primaryHeadroom
    targetLambda
    options.ExtraAmbientSpd = 0;
end

%% Preparations for the calculations.
%
% This part prepares variables for following calculations.
%
% Get background Spd and LMS. Pass this through quantization
bgPrimariesQuantized = SettingsToPrimary(subprimaryCalstructData,PrimaryToSettings(subprimaryCalstructData,bgPrimaries));
bgSpd = PrimaryToSpd(subprimaryCalstructData,bgPrimariesQuantized) + ...
        options.ExtraAmbientSpd;
bgLMS = T_cones * bgSpd;

% Set the target contrast maximum value with headroom.
targetContrastReMaxWithHeadroom = targetPrimaryHeadroom * targetContrastReMax;

%% Get primaries based on contrast specification.
targetLMSContrast = targetContrastReMaxWithHeadroom * targetMaxLMSContrast;
[isolatingModulationPrimaries] = ReceptorIsolateSpectral(T_cones,targetLMSContrast,subprimaryCalstructData.get('P_device'),bgPrimaries,bgPrimaries, ...
    primaryHeadroom,B_natural,projectIndices,targetLambda,subprimaryCalstructData.get('P_ambient')+options.ExtraAmbientSpd,'EXCITATIONS',false);
isolatingPrimaries = isolatingModulationPrimaries + bgPrimaries;

% Quantize.
isolatingPrimariesQuantized = SettingsToPrimary(subprimaryCalstructData,PrimaryToSettings(subprimaryCalstructData,isolatingPrimaries));

% Get isolating spd and LMS contrasts.
isolatingSpd = PrimaryToSpd(subprimaryCalstructData,isolatingPrimariesQuantized);
isolatingSpdForContrastCalc = isolatingSpd + options.ExtraAmbientSpd;
isolatingLMSForContrastCalc = T_cones * isolatingSpdForContrastCalc;
isolatingContrast = ExcitationsToContrast(isolatingLMSForContrastCalc,bgLMS);

% Report.
fprintf('Desired/obtained contrasts:\n');
for rr = 1:length(targetLMSContrast)
    fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,targetLMSContrast(rr),isolatingContrast(rr));
end
fprintf('Min/max primaries: %0.4f, %0.4f\n',min(isolatingPrimariesQuantized),max(isolatingPrimariesQuantized));

end