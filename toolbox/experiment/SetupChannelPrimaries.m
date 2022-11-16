function [screenPrimaryChannelObject,backgroundChannelObject] = SetupChannelPrimaries(colorDirectionParams,channelCalObjs,projectIndices,options)
% Set up channel primaries that reproduce desired cone contrasts.
%
% Syntax:
%    [screenPrimaryChannelObj,bgChannelObject] = SetupChannelPrimaries(colorDirectionParams,channelCalObjs,projectIndices)
%
% Description:
%    This finds the channel primareis that reproduce desired cone
%    contrasts. Here we find the primaries of the background first that
%    matches a specific chromaticity, then we find the channel primaries
%    after.
%
% Inputs:
%    colorDirectionParams          - Structure with the parameters to
%                                    calculate a contrast gabor image.
%    channelCalObjs                - Channel calibration objects.
%    projectIndices                - Wavelength range that is used to
%                                    enforce the smoothness through the
%                                    projection onto an underlying basis
%                                    set.
%
% Outputs:
%    screenPrimaryChannelObject    - Structure that contains the results of
%                                    the target screen primary channel
%                                    primaries.
%    backgroundChannelObject       - Structure that contains the background
%                                    channel primaries.
%
% Optional key/value pairs:
%    verbose                       - Boolean. Default true. Controls
%                                    plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,ga,smo           - Wrote it.
%   01/24/22  smo                  - Made it work.
%   11/16/22  smo                  - Update the target lambda right.

%% Set parameters.
arguments
    colorDirectionParams
    channelCalObjs
    projectIndices
    options.verbose (1,1) = true
end

%% Find background primaries to acheive desired xy at intensity scale of display.
%
% Set parameters for getting desired background primaries.
bgPrimaryHeadRoom = 0;
bgTargetLambda = colorDirectionParams.targetLambda;
targetBgXYZ = xyYToXYZ([colorDirectionParams.targetBgxy ; 1]);
nScreenPrimaries = size(colorDirectionParams.channelCalNames,2);

% Make a loop for getting background for all primaries.
%
% Passing true for key 'Scale' causes these to be scaled reasonably
% relative to gamut, which is why we can set the target luminance
% arbitrarily to 1 just above. The scale factor determines where in the
% approximate channel gamut we aim the background at.
for pp = 1:nScreenPrimaries
    [channelBackgroundPrimaries(:,pp),channelBackgroundSpd(:,pp),channelBackgroundXYZ(:,pp)] = ...
        FindBgChannelPrimaries(targetBgXYZ, colorDirectionParams.T_xyz, channelCalObjs{pp}, colorDirectionParams.B_natural{pp}, ...
        projectIndices, bgPrimaryHeadRoom, bgTargetLambda, 'scaleFactor', 0.6, 'Scale', true, 'Verbose', options.verbose);
end

% Check if all primaries are within the range from 0 to 1.
if (any(channelBackgroundPrimaries < 0) | any(channelBackgroundPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end
if (options.verbose)
fprintf('Background primary min: %0.2f, max: %0.2f, mean: %0.2f\n', ...
    min(channelBackgroundPrimaries(:)), max(channelBackgroundPrimaries(:)), mean(channelBackgroundPrimaries(:)));
end

%% Get isolating primaries for all screen primaries.
for pp = 1:nScreenPrimaries
    % The ambient with respect to which we compute contrast is from all
    % three primaries, which we handle via the extraAmbientSpd key-value
    % pair in the call.  The extra is for the primaries not being found in
    % the current call - the contribution from the current primary is known
    % because we pass the primaries for the background.
    otherPrimaries = setdiff(1:nScreenPrimaries,pp);
    extraAmbientSpd = 0;
    for oo = 1:length(otherPrimaries)
        extraAmbientSpd = extraAmbientSpd + channelBackgroundSpd(:,otherPrimaries(oo));
    end
    
    % Get isolating screen primaries.
    [screenPrimaryPrimaries(:,pp),screenPrimaryPrimariesQuantized(:,pp),screenPrimarySpd(:,pp),screenPrimaryContrast(:,pp),screenPrimaryModulationPrimaries(:,pp)] ...
        = FindChannelPrimaries(colorDirectionParams.targetScreenPrimaryContrastDir(:,pp), ...
        colorDirectionParams.targetPrimaryHeadroom, colorDirectionParams.targetScreenPrimaryContrasts(pp), channelBackgroundPrimaries(:,pp), ...
        colorDirectionParams.T_cones, channelCalObjs{pp}, colorDirectionParams.B_natural{pp}, projectIndices, colorDirectionParams.primaryHeadroom,...
        colorDirectionParams.targetLambda, 'ExtraAmbientSpd', extraAmbientSpd);
    
    % We can wonder about how close to gamut our primaries are.  Compute
    % that here.
    primaryGamutScaleFactor(pp) = MaximizeGamutContrast(screenPrimaryModulationPrimaries(:,pp),channelBackgroundPrimaries(:,pp));
    if(options.verbose)
        fprintf('\tPrimary %d, gamut scale factor is %0.3f\n',pp,primaryGamutScaleFactor(pp));
    end
    
    % Find the channel settings that correspond to the desired screen
    % primaries.
    screenPrimarySettings(:,pp) = PrimaryToSettings(channelCalObjs{pp},screenPrimaryPrimaries(:,pp));
end

%% How close are spectra to subspace defined by basis?
for pp = 1:nScreenPrimaries
    isolatingNaturalApproxSpd(:,pp) = colorDirectionParams.B_natural{pp} * (colorDirectionParams.B_natural{pp}(projectIndices,:)\screenPrimarySpd(projectIndices,pp));
end

% Plot it.
if (options.verbose)
    figure; clf;
    for pp = 1:nScreenPrimaries
        subplot(2,2,pp); hold on;
        plot(colorDirectionParams.wls,screenPrimarySpd(:,pp),'b','LineWidth',2);
        plot(colorDirectionParams.wls,isolatingNaturalApproxSpd(:,pp),'r:','LineWidth',1);
        plot(colorDirectionParams.wls(projectIndices),screenPrimarySpd(projectIndices,pp),'b','LineWidth',4);
        plot(colorDirectionParams.wls(projectIndices),isolatingNaturalApproxSpd(projectIndices,pp),'r:','LineWidth',3);
        xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
        title(append('Primary ', num2str(pp)));
    end
end

%% Save the results in the structures.
%
% Screen background priamry struct.
backgroundChannelObject.channelBackgroundPrimaries = channelBackgroundPrimaries;
backgroundChannelObject.channelBackgroundSpd = channelBackgroundSpd;
backgroundChannelObject.channelBackgroundXYZ = channelBackgroundXYZ;

% Screen primary channel struct.
screenPrimaryChannelObject.screenPrimaryPrimaries= screenPrimaryPrimaries;
screenPrimaryChannelObject.screenPrimaryPrimariesQuantized = screenPrimaryPrimariesQuantized;
screenPrimaryChannelObject.screenPrimarySpd = screenPrimarySpd;
screenPrimaryChannelObject.screenPrimaryContrast = screenPrimaryContrast;
screenPrimaryChannelObject.screenPrimaryModulationPrimaries = screenPrimaryModulationPrimaries;
screenPrimaryChannelObject.screenPrimarySettings = screenPrimarySettings;
screenPrimaryChannelObject.isolatingNaturalApproxSpd = isolatingNaturalApproxSpd;

end
