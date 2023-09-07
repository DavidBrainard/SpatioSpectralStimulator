function [backgroundScreenPrimaryObject] = SetupBackground(colorDirectionParams,screenCalObj,backgroundChannelObject,options)
% Find the background screen primaries that reproduces a desired chromaticity.
%
% Syntax:
%    [backgroundScreenPrimaryObject] = SetupBackground(colorDirectionParams,screenCalObj,backgroundChannelObject)
%
% Description:
%    This calculates the background screen primaries that reproduces a
%    desired chromaticity. As we decided the channel primaries, here we
%    find the screen primaries for the background that reproduces the same
%    chromaticity as the one we found when calculating the channel
%    primaries.
%
% Inputs:
%    colorDirectionParams              - Structure with the parameters to
%                                        calculate a contrast gabor image.
%    screenCalObj                      - Screen calibration object.
%    backgroundChannelObject           - Structure that contains the
%                                        background channel primaries.
%
% Outputs:
%    backgroundScreenPrimaryObject     - Structure with the background
%                                        screen primaries that reproduces a
%                                        specific chromaticity.
%
% Optional key/value pairs:
%    screenBackgroundScaleFactor       - Adjust this to keep the background
%                                        in gamut.
%    verbose                           - Boolean. Default true. Controls
%                                        plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,gka,smo              - Wrote it.
%   01/24/22  smo                      - Made it work.

%% Set parameters.
arguments
    colorDirectionParams
    screenCalObj
    backgroundChannelObject
    options.screenBackgroundScaleFactor (1,1) = 0.5
    options.verbose (1,1) = true
end

%% Set up desired background.
%
% We aim for the background that we said we wanted when we built the screen primaries.
if (isfield(colorDirectionParams,'T_receptors'))
     desiredBgExcitations = options.screenBackgroundScaleFactor * colorDirectionParams.T_receptors * sum(backgroundChannelObject.channelBackgroundSpd,2);
else
    desiredBgExcitations = options.screenBackgroundScaleFactor * colorDirectionParams.T_cones * sum(backgroundChannelObject.channelBackgroundSpd,2);
end

screenBgSettings = SensorToSettings(screenCalObj,desiredBgExcitations);
screenBgExcitations = SettingsToSensor(screenCalObj,screenBgSettings);

% Plot it.
if (options.verbose)
    figure; clf; hold on;
    plot(desiredBgExcitations,screenBgExcitations,'ro','MarkerFaceColor','r','MarkerSize',12);
    axis('square');
    xlim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
    ylim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
    xlabel('Desired bg excitations'); ylabel('Obtained bg excitations');
    title('Check that we obtrain desired background excitations');
    fprintf('Screen settings to obtain background: %0.2f, %0.2f, %0.2f\n', ...
        screenBgSettings(1),screenBgSettings(2),screenBgSettings(3));
end

%% Save the results in a struct.
backgroundScreenPrimaryObject.desiredBgExcitations = desiredBgExcitations;
backgroundScreenPrimaryObject.screenBgSettings = screenBgSettings;
backgroundScreenPrimaryObject.screenBgExcitations = screenBgExcitations;

end
