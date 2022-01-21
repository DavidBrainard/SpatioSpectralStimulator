function [] = SetupChannelPrimaries()
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

%% Find background primaries to acheive desired xy at intensity scale of display.
%
% Set parameters for getting desired background primaries.
primaryHeadRoom = 0;
targetLambda = 3;
targetBgXYZ = xyYToXYZ([colorDirectionParams.targetBgxy ; 1]);

% Adjust these to keep background in gamut
% primaryBackgroundScaleFactor = 0.5;
screenBackgroundScaleFactor = 0.5;

% Make a loop for getting background for all primaries.
% Passing true for key 'Scale' causes these to be scaled reasonably
% relative to gamut, which is why we can set the target luminance
% arbitrarily to 1 just above. The scale factor determines where in the
% approximate channel gamut we aim the background at.
for pp = 1:nScreenPrimaries
    [channelBackgroundPrimaries(:,pp),channelBackgroundSpd(:,pp),channelBackgroundXYZ(:,pp)] = ...
        FindBgChannelPrimaries(targetBgXYZ, colorDirectionParams.T_xyz, channelCalObjs{pp}, colorDirectionParams.B_natural{pp}, ...
        projectIndices, primaryHeadRoom, targetLambda, 'scaleFactor', 0.6, 'Scale', true, 'Verbose', VERBOSE);
end
if (any(channelBackgroundPrimaries < 0) | any(channelBackgroundPrimaries > 1))
    error('Oops - primaries should always be between 0 and 1');
end
fprintf('Background primary min: %0.2f, max: %0.2f, mean: %0.2f\n', ...
    min(channelBackgroundPrimaries(:)), max(channelBackgroundPrimaries(:)), mean(channelBackgroundPrimaries(:)));

% Get isolating primaries for all screen primaries.
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
        colorDirectionParams.targetPrimaryHeadroom,colorDirectionParams.targetScreenPrimaryContrasts(pp),channelBackgroundPrimaries(:,pp), ...
        colorDirectionParams.T_cones,channelCalObjs{pp},colorDirectionParams.B_natural{pp},projectIndices,colorDirectionParams.primaryHeadroom,...
        colorDirectionParams.targetLambda,'ExtraAmbientSpd',extraAmbientSpd);
    
    % We can wonder about how close to gamut our primaries are.  Compute
    % that here.
    primaryGamutScaleFactor(pp) = MaximizeGamutContrast(screenPrimaryModulationPrimaries(:,pp),channelBackgroundPrimaries(:,pp));
    fprintf('\tPrimary %d, gamut scale factor is %0.3f\n',pp,primaryGamutScaleFactor(pp));
    
    % Find the channel settings that correspond to the desired screen
    % primaries.
    screenPrimarySettings(:,pp) = PrimaryToSettings(channelCalObjs{pp},screenPrimaryPrimaries(:,pp));
end

%% How close are spectra to subspace defined by basis?
%
% This part has been updated using the loop to make it short.
for pp = 1:nScreenPrimaries
    isolatingNaturalApproxSpd(:,pp) = colorDirectionParams.B_natural{pp} * (colorDirectionParams.B_natural{pp}(projectIndices,:)\screenPrimarySpd(projectIndices,pp));
end

% Plot of the screen primary spectra.
figure; clf;
for pp = 1:nScreenPrimaries
    subplot(2,2,pp); hold on;
    plot(wls,screenPrimarySpd(:,pp),'b','LineWidth',2);
    plot(wls,isolatingNaturalApproxSpd(:,pp),'r:','LineWidth',1);
    plot(wls(projectIndices),screenPrimarySpd(projectIndices,pp),'b','LineWidth',4);
    plot(wls(projectIndices),isolatingNaturalApproxSpd(projectIndices,pp),'r:','LineWidth',3);
    xlabel('Wavelength (nm)'); ylabel('Power (arb units)');
    title(append('Primary ', num2str(pp)));
end
end
