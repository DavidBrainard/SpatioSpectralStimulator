function [isolatingSpdSet,isolatingContrastSet] = FindDesiredTargetPrimaries(T_cones,B_natural,subprimaryCalstructData, ...
                                                  targetMaxLMSContrast,targetPrimaryHeadroom,targetContrastReMax,bgPrimaries,projectIndices,primaryHeadroom,targetLambda)

% Find the desired target primaries so that the three primaries can produce
% desired LMS contrast.
%
% Syntax:
%
%
% Description:
%    This function searches for three primaries based on the desired LMS
%    contrast from the background. The shape of the final obtained spectrums will
%    depend on the variable 'targetLambda' which constrains the smoothness
%    of the curve.
%
% Inputs:
%    T_cones -                    Spectral cone sensitivities in standard PTB format.
%    B_natural -                  Basis functions to enforce the smoothness
%                                 on the desired background spectrum.
%                                 This also should be in standard PTB
%                                 format.
%    subPrimaryCalstructData -    The subprimary calstructdata to calculate
%                                 some routines including 'PrimaryToSpd','SettingsToPrimary',
%                                 and 'PrimaryToSettings'
%    targetMaxLMSContrast -       The maximum LMS contrast values for the target 
%                                 are stored here for all three primaries.
%                                 Each column should be aligned as each
%                                 primary. 
%    targetPrimaryHeadroom -      
%    targetContrastReMax -        
%    bgPrimaries -                
%    projectIndices -             Define wavelength range that will be used to enforce the smoothnes
%                                 thorugh the projection onto an underlying basis set.
%    primaryHeadRoom -            Set the headroom to prevent running into
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
%    isolatingSpdSet -            This contains the final spectrum results 
%                                 of the three primaries in struct
%                                 variable. 
%                                 (cf. isolatingSpdSet = {primary1Spd primary2Spd primary3Spd}).
%    isolatingContrastSet -       This contains the final obtained contrast
%                                 results for three primaries.
%                                 (cf. isolatingContrastSet = {primary1Contrast primary2Contrast primary3Contrast}).
%
% Optional key/value pairs:
%    N/A
%
% History:
%    09/29/21  smo  Started on it

%% Preparations for the calculations.
%
% This part prepares variables for following calculations.
%
% Get background Spd and LMS.
primarySettingForBackground = 1; % Choose which primary to use for the backgroud. [1,2,3]. / For now, it is set as 1 arbitrarily.
bgSpd = PrimaryToSpd(subprimaryCalstructData{primarySettingForBackground},bgPrimaries); % This should be updated later for using three primaries for the background?
bgLMS = T_cones*bgSpd;

% Set the target contrast maximum value with headroom.
targetContrastReMaxWithHeadroom = targetPrimaryHeadroom*targetContrastReMax;

%% Get primaries based on contrast specification.
%
% Make a loop for three primaries.
numPrimaries = 3;
for pp = 1:numPrimaries;
    % Set the target contrast per each primary.
    targetLMSContrast_temp = targetContrastReMaxWithHeadroom*targetMaxLMSContrast(:,pp);
    [isolatingModulationPrimaries_temp] = ReceptorIsolateSpectral(T_cones,targetLMSContrast_temp,subprimaryCalstructData{pp}.get('P_device'),bgPrimaries,bgPrimaries, ...
                                          primaryHeadRoom,B_natural,projectIndices,targetLambda,subprimaryCalstructData{pp}.get('P_ambient'),'EXCITATIONS',false);
    isolatingPrimaries_temp = isolatingModulationPrimaries_temp + bgPrimaries;
    
    % Quantize.
    isolatingPrimariesQuantized_temp = SettingsToPrimary(subprimaryCalstructData{pp},PrimaryToSettings(subprimaryCalstructData{pp},isolatingPrimaries_temp));
    
    % Report.
    isolatingSpd_temp = PrimaryToSpd(subprimaryCalstructData{pp},isolatingPrimariesQuantized_temp);
    isolatingLMS_temp = T_cones*isolatingSpd_temp;
    isolatingContrast_temp = ExcitationsToContrast(isolatingLMS_temp,bgLMS);
    fprintf('Desired/obtained contrasts 1\n');
    for rr = 1:length(targetLMSContrast_temp)
        fprintf('\tReceptor %d (desired/obtained): %0.3f, %0.3f\n',rr,targetLMSContrast_temp(rr),isolatingContrast_temp(rr));
    end
    fprintf('Min/max primaries 1: %0.4f, %0.4f\n',min(isolatingPrimariesQuantized_temp),max(isolatingPrimariesQuantized_temp));
    
    % Save the target isolating Spd and Contrast.
    variableName_temp = append('isolatingSpd',num2str(pp));
    isolatingSpdSet.variableName_temp = isolatingSpd_temp;
    isolatingContrastSet{pp} = isolatingContrast_temp;
end

end