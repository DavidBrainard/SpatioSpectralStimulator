function [ISETBioGaborObject] = MakeISETBioSceneFromImage(colorDirectionParams,gaborImageObject,standardGaborCalObject,...
    ISETBioDisplayObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg,options)
% Make ISETBio scene from the gabor image.
%
% Syntax:
%    [ISETBioGaborCalObject] = MakeISETBioSceneFromImage(colorDirectionParams,gaborImageObject,standardGaborCalObject,...
%                              ISETBioDisplayObject,stimulusHorizSizeMeters,stimulusHorizSizeDeg)
%
% Description:
%    This puts the target gabor image into ISETBio scene. 
%
% Inputs:
%    colorDirectionParams          - Structure with the parameters to
%                                    calculate a contrast gabor image.
%    gaborImageObject              - Structure with the gabor contrast image in
%                                    image format.
%    standardGaborCalObject        - Structure with the gabor contrasts
%                                    and settings in a cal format.
%    ISETBioDisplayObject          - Structure with the parameters to make the
%                                    ISETBio scene from image.
%    stimulusHorizSizeMeters       - The horizontal size of the gabor image
%                                    in meters.
%    stimulusHorizSizeDeg          - The horizontal size of the gabor image
%                                    in degrees.
%
% Outputs:
%    ISETBioGaborObject            - Structure with gabor contrast image in
%                                    the ISETBio scene format.
%
% Optional key/value pairs:
%    verbose                       - Boolean. Default true. Controls
%                                    plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio, GetSettingsFromISETBioScene

% History:
%   01/21/22  dhb,gka,smo     - Wrote it.
%   01/24/22  smo             - Made it work.

%% Set parameters.
arguments
    colorDirectionParams
    gaborImageObject
    standardGaborCalObject
    ISETBioDisplayObject
    stimulusHorizSizeMeters
    stimulusHorizSizeDeg
    options.verbose (1,1) = true
end

%% Put the image into an ISETBio scene.
%
% These calls are a bit slow for large images and the fine wavelength
% sampling used here. But these would be done as pre-compute steps so
% it doesn't seem worth trying to optimize at this point.
ISETBioGaborScene = sceneFromFile(gaborImageObject.standardSettingsGaborImage,'rgb', [], ISETBioDisplayObject);

% Show the image on ISETBio scene window.
if (options.verbose)
    sceneWindow(ISETBioGaborScene);
end

% Check stimulus dimensions match. These are good to about a percent, which
% we can live with.
stimulusHorizSizeMetersChk = sceneGet(ISETBioGaborScene,'width');
stimulusHorizSizeDegChk = sceneGet(ISETBioGaborScene,'horizontal fov');
if (abs(stimulusHorizSizeMeters - stimulusHorizSizeMetersChk)/stimulusHorizSizeMeters > 0.01)
    error('Horizontal size in meters mismatch of too much');
end
if (abs(stimulusHorizSizeDeg - stimulusHorizSizeDegChk)/stimulusHorizSizeDeg > 0.01)
    error('Horizontal size in deg mismatch of too much');
end

%% Calculate cone excitations from the ISETBio scene.
% These should match what we get when we compute
% outside of ISETBio. And indeed!
%
% ISETBio energy comes back as power per nm, we need to convert to power
% per wlband to work with PTB, by multiplying by S(2).
ISETBioGaborImage = sceneGet(ISETBioGaborScene,'energy') * colorDirectionParams.S(2);
[ISETBioGaborCal,ISETBioM,ISETBioN] = ImageToCalFormat(ISETBioGaborImage);
ISETBioPredictedExcitationsGaborCal = colorDirectionParams.T_cones * ISETBioGaborCal;
limMin = 0.01; limMax = 0.02;

% Plot it to comapare the cone excitations between before and after passing
% the ISETBio scene.
if (options.verbose)
    figure; clf; hold on;
    plot(standardGaborCalObject.standardPredictedExcitationsGaborCal(1,:), ISETBioPredictedExcitationsGaborCal(1,:),'r+');
    plot(standardGaborCalObject.standardPredictedExcitationsGaborCal(2,:), ISETBioPredictedExcitationsGaborCal(2,:),'g+');
    plot(standardGaborCalObject.standardPredictedExcitationsGaborCal(3,:), ISETBioPredictedExcitationsGaborCal(3,:),'b+');
    plot([limMin limMax], [limMin limMax]);
    xlabel('Standard Cone Excitations');
    ylabel('ISETBio Cone Excitations');
    axis('square'); xlim([limMin limMax]); ylim([limMin limMax]);
    title('Cone Excitations Comparison');
end

% Check if it predicts well.
if (max(abs(standardGaborCalObject.standardPredictedExcitationsGaborCal(:) - ISETBioPredictedExcitationsGaborCal(:)) ./ ...
        standardGaborCalObject.standardPredictedExcitationsGaborCal(:)) > 1e-6)
    error('Standard and ISETBio data do not agree well enough');
end

% Print out if everything goes well.
if (options.verbose)
    disp('Gabor image has been successfully calculated from the ISETBio scene!');
end

%% Save the results in a struct.
ISETBioGaborObject.ISETBioGaborScene = ISETBioGaborScene;
ISETBioGaborObject.ISETBioGaborImage = ISETBioGaborImage;
ISETBioGaborObject.ISETBioPredictedExcitationsGaborCal = ISETBioPredictedExcitationsGaborCal;
ISETBioGaborObject.ISETBioGaborCal = ISETBioGaborCal;

end
