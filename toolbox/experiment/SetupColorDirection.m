function colorDirectionParams = SetupColorDirection(conditionName,options)
% Set up a parameter structure for a color direction.
%
% Syntax:
%    colorDirectionParams = SetupColorDirection(conditionName)
%
% Description:
%    Take in a condition name string and return a structure with the
%    parameters required to compute stimuli for that condition. This then
%    allows computation of stimuli to present.
%
%    How these parameters get used is illustrated in tutorial form
%    elsewhere, see the "See also" list below.
%
% Inputs:
%    conditionName                - String containing the type name for the
%                                   direction.
%
% Outputs:
%    colorDirectionParams         - Structure with the needed parameters.
%
% Optional key/value pairs:
%    spatialGaborTargetContrast   - Maximum contrast in the target
%                                   gabor image that we will make.
%    targetLambda                 - Default to 3. Controls the smootheness
%                                   of the target primaries. Higher value
%                                   means more constrained to make it
%                                   smoother.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/18/22  dhb, smo             - Wrote it
%   01/26/22  smo                  - Added an option to set the target
%                                    contrast of the gabor image.
%   11/15/22  smo                  - Added an option to set target lambda.

%% Set parameters.
arguments
    conditionName {mustBeMember(conditionName,{'LminusMSmooth','ConeIsolating'})}
    options.spatialGaborTargetContrast (1,1) = 0.04
    options.targetScreenPrimaryContrasts (1,1) = 0.05
    options.targetLambda (1,1) = 3
end

%% Set some initial parameters on the struct here.
%
% Some parameters which might become key/value pairs.
colorDirectionParams.fieldSizeDegrees = 2;

% Set condition name.
colorDirectionParams.conditionName = conditionName;

% Set wavelength support.
colorDirectionParams.S = [380 2 201];
colorDirectionParams.wls = SToWls(colorDirectionParams.S);

%% Set key stimulus parameters
%
% Condition Name.
switch (colorDirectionParams.conditionName)
    case 'LminusMSmooth'
        % Background xy.
        %
        % Specify the chromaticity, but we'll chose the luminance based
        % on the range available in the device.
        colorDirectionParams.targetBgxy = [0.3127 0.3290]';
        
        % Target color direction and max contrasts.
        %
        % This is the basic desired modulation direction positive excursion. We go
        % equally in positive and negative directions.  Make this unit vector
        % length, as that is good convention for contrast.
        colorDirectionParams.targetStimulusContrastDir = [1 -1 0]';
        colorDirectionParams.targetStimulusContrastDir = colorDirectionParams.targetStimulusContrastDir / norm(colorDirectionParams.targetStimulusContrastDir);
        
        % Specify desired primary properties.
        %
        % These are the target contrasts for the three primaries. We want these to
        % span a triangle around the line specified above. Here we define that
        % triangle by hand.  May need a little fussing for other directions, and
        % might be able to autocompute good choices.
        colorDirectionParams.targetScreenPrimaryContrastDir(:,1) = [-1 1 0]';
        colorDirectionParams.targetScreenPrimaryContrastDir(:,1) = colorDirectionParams.targetScreenPrimaryContrastDir(:,1) / norm(colorDirectionParams.targetScreenPrimaryContrastDir(:,1));
        colorDirectionParams.targetScreenPrimaryContrastDir(:,2) = [1 -1 0.5]';
        colorDirectionParams.targetScreenPrimaryContrastDir(:,2) = colorDirectionParams.targetScreenPrimaryContrastDir(:,2) / norm(colorDirectionParams.targetScreenPrimaryContrastDir(:,2));
        colorDirectionParams.targetScreenPrimaryContrastDir(:,3) = [1 -1 -0.5]';
        colorDirectionParams.targetScreenPrimaryContrastDir(:,3) = colorDirectionParams.targetScreenPrimaryContrastDir(:,3) / norm(colorDirectionParams.targetScreenPrimaryContrastDir(:,3));
        
        % Set parameters for getting desired target primaries.
        colorDirectionParams.targetScreenPrimaryContrasts = ones(1,3) * options.targetScreenPrimaryContrasts;
        colorDirectionParams.targetPrimaryHeadroom = 1.05;
        colorDirectionParams.primaryHeadroom = 0;
        colorDirectionParams.targetLambda = options.targetLambda;
        
        % We may not need the whole direction contrast excursion. Specify max
        % contrast we want relative to that direction vector.
        % The first number is
        % the amount we want to use, the second has a little headroom so we don't
        % run into numerical error at the edges. The second number is used when
        % defining the three primaries, the first when computing desired weights on
        % the primaries.
        colorDirectionParams.spatialGaborTargetContrast = options.spatialGaborTargetContrast;
        
        % Set up basis to try to keep spectra close to.
        %
        % This is how we enforce a smoothness or other constraint
        % on the spectra.  What happens in the routine that finds
        % primaries is that there is a weighted error term that tries to
        % maximize the projection onto a passed basis set.
        colorDirectionParams.basisType = 'fourier';
        colorDirectionParams.nFourierBases = 7;
        switch (colorDirectionParams.basisType)
            case 'cieday'
                load B_cieday
                colorDirectionParams.B_naturalRaw = SplineSpd(S_cieday, B_cieday, colorDirectionParams.S);
            case 'fourier'
                colorDirectionParams.B_naturalRaw = MakeFourierBasis(colorDirectionParams.S, colorDirectionParams.nFourierBases);
            otherwise
                error('Unknown basis set specified');
        end
        colorDirectionParams.B_natural{1} = colorDirectionParams.B_naturalRaw;
        colorDirectionParams.B_natural{2} = colorDirectionParams.B_naturalRaw;
        colorDirectionParams.B_natural{3} = colorDirectionParams.B_naturalRaw;
        
    otherwise
        error('Unknown condition name specified')
end

%% Direction independent device parameters.
%
% Calibration filenames/params.
%
% This is a standard calibration file for the DLP screen,
% with the channels set to something.
colorDirectionParams.screenCalName = 'SACC';
colorDirectionParams.screenNInputLevels = 256;

% These are the calibration files for each of the screen primaries, which
% then entails measuring the spectra of all the channels for that
% screen primary.
colorDirectionParams.channelCalNames = {'SACCPrimary1' 'SACCPrimary2' 'SACCPrimary3'};
colorDirectionParams.channelNInputLevels = 253;

%% Direction independent colorimetery parameters.
%
% Cone fundamentals and XYZ CMFs.
colorDirectionParams.psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
colorDirectionParams.psiParamsStruct.coneParams.fieldSizeDegrees = colorDirectionParams.fieldSizeDegrees;
colorDirectionParams.T_cones = ComputeObserverFundamentals(colorDirectionParams.psiParamsStruct.coneParams,colorDirectionParams.S);

% Judd-Vos XYZ color matching function.
load T_xyzJuddVos;
colorDirectionParams.T_xyz = SplineCmf(S_xyzJuddVos, 683*T_xyzJuddVos, colorDirectionParams.S);

end
