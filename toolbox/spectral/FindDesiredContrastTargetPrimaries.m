function [targetSpdSet,targetContrastSet] = FindDesiredContrastTargetPrimaries(targetMaxLMSContrast,targetPrimaryHeadroom,targetContrastReMax,bgPrimaries, ...
    T_cones,subprimaryCalstructData,B_natural,projectIndices,primaryHeadroom,targetLambda)

% Find the desired target primaries so that the three primaries can produce
% desired LMS contrast.
%
% Syntax:
%    [targetSpdSet,targetContrastSet] = FindDesiredContrastTargetPrimaries(targetMaxLMSContrast,targetPrimaryHeadroom,targetContrastReMax,bgPrimaries, ...
%                                       T_cones,subprimaryCalstructData,B_natural,projectIndices,primaryHeadroom,targetLambda)
%
% Description:
%    This function searches for three primaries based on the desired LMS
%    contrast from the background. The shape of the final obtained spectrums will
%    depend on the variable 'targetLambda' which constrains the smoothness
%    of the curve.
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
%    targetSpdSet -               This contains the final spectrum results
%                                 of the three primaries in struct
%                                 variable.
%                                 (cf. targetSpdSet = {primary1Spd primary2Spd primary3Spd}).
%    targetContrastSet -          This contains the final obtained contrast
%                                 results for three primaries.
%                                 (cf. targetContrastSet = {primary1Contrast primary2Contrast primary3Contrast}).
%
% Optional key/value pairs:
%    N/A
%
% History:
%    09/29/21  smo  Started on it

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
end

%% Preparations for the calculations.
%
% This part prepares variables for following calculations.
%
% Get background Spd and LMS.
primarySettingForBackground = 1; % Choose which primary to use for the backgroud. / For now, it is set as 1 arbitrarily.
bgSpd = PrimaryToSpd(subprimaryCalstructData{primarySettingForBackground},bgPrimaries); % This should be updated later for using three primaries for the background?
bgLMS = T_cones * bgSpd;

% Set the target contrast maximum value with headroom.
targetContrastReMaxWithHeadroom = targetPrimaryHeadroom * targetContrastReMax;

%% Get primaries based on contrast specification.
%
% Make a loop for three primaries.
numPrimaries = length(subprimaryCalstructData);
for pp = 1:numPrimaries;
    % Set the target contrast per each primary.
    targetLMSContrast_temp = targetContrastReMaxWithHeadroom * targetMaxLMSContrast(:,pp);
    [isolatingModulationPrimaries_temp] = ReceptorIsolateSpectral(T_cones,targetLMSContrast_temp,subprimaryCalstructData{pp}.get('P_device'),bgPrimaries,bgPrimaries, ...
        primaryHeadroom,B_natural,projectIndices,targetLambda,subprimaryCalstructData{pp}.get('P_ambient'),'EXCITATIONS',false);
    isolatingPrimaries_temp = isolatingModulationPrimaries_temp + bgPrimaries;
    
    % Quantize.
    isolatingPrimariesQuantized_temp = SettingsToPrimary(subprimaryCalstructData{pp},PrimaryToSettings(subprimaryCalstructData{pp},isolatingPrimaries_temp));
    
    % Report.
    isolatingSpd_temp = PrimaryToSpd(subprimaryCalstructData{pp},isolatingPrimariesQuantized_temp);
    isolatingLMS_temp = T_cones * isolatingSpd_temp;
    isolatingContrast_temp = ExcitationsToContrast(isolatingLMS_temp,bgLMS);
    fprintf('Desired/obtained contrasts 1\n');
    for rr = 1:length(targetLMSContrast_temp)
        fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,targetLMSContrast_temp(rr),isolatingContrast_temp(rr));
    end
    fprintf('Min/max primaries 1: %0.4f, %0.4f\n',min(isolatingPrimariesQuantized_temp),max(isolatingPrimariesQuantized_temp));
    
    % Save the target isolating Spd and Contrast.
    % It saves the data temporarily and the data will be stored in the
    % 'struct' in the following lines after the loop.
    isolatingSpdSet_temp(:,pp) = isolatingSpd_temp;
    isolatingContrastSet_temp(:,pp) = isolatingContrast_temp;
end

% Make a struct for each result of the Spd and the Contrast. There
% should be better way to do it, so this part would be changed later.
targetSpdSet = struct('isolatingSpd1',isolatingSpdSet_temp(:,1),'isolatingSpd2',isolatingSpdSet_temp(:,2),'isolatingSpd3',isolatingSpdSet_temp(:,3));
targetContrastSet = struct('isolatingContrast1',isolatingContrastSet_temp(:,1),'isolatingContrast2',isolatingContrastSet_temp(:,2),'isolatingContrast3',isolatingContrastSet_temp(:,3));
end