function [bgScreenPrimaryObj] = SetupBackground(colorDirectionParams,screenCalObj,bgChannelObject)
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
    colorDirectionParams
    screenCalObj
    bgChannelObject
end

%% Set up desired background.
%
% Adjust these to keep background in gamut
screenBackgroundScaleFactor = 0.5;

% We aim for the background that we said we wanted when we built the screen primaries.
desiredBgExcitations = screenBackgroundScaleFactor * colorDirectionParams.T_cones * sum(bgChannelObject.channelBackgroundSpd,2);
screenBgSettings = SensorToSettings(screenCalObj,desiredBgExcitations);
screenBgExcitations = SettingsToSensor(screenCalObj,screenBgSettings);

% Plot it.
figure; clf; hold on;
plot(desiredBgExcitations,screenBgExcitations,'ro','MarkerFaceColor','r','MarkerSize',12);
axis('square');
xlim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
ylim([min([desiredBgExcitations ; screenBgExcitations]),max([desiredBgExcitations ; screenBgExcitations])]);
xlabel('Desired bg excitations'); ylabel('Obtained bg excitations');
title('Check that we obtrain desired background excitations');
fprintf('Screen settings to obtain background: %0.2f, %0.2f, %0.2f\n', ...
    screenBgSettings(1),screenBgSettings(2),screenBgSettings(3));

%% Save the results in a struct.
bgScreenPrimaryObj.desiredBgExcitations = desiredBgExcitations;
bgScreenPrimaryObj.screenBgSettings = screenBgSettings;
bgScreenPrimaryObj.screenBgExcitations = screenBgExcitations;

end

