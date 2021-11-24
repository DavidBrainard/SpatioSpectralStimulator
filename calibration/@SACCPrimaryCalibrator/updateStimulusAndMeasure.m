% Method to update the stimulus and conduct a single radiometric measurement.
%
% History:
%   11/24/2021 smo  Now measurement uses SACC measurement function instead
%                   of radiometerOBJ.

function [measurement, S] = updateStimulusAndMeasure(obj, bgSettings, targetSettings, useBitsPP)

% Display current subprimary settings.
if (obj.options.verbosity > 1)
    for i=1:obj.nSubprimaries
        fprintf('        Target settings %2.0f: %2.0f \n',i,round((obj.nInputLevels-1)*targetSettings(i)));
    end
end

% Update background and target stimuli.
obj.updateBackgroundAndTarget(bgSettings, targetSettings, useBitsPP);

% Make a delay before the measurement for warming-up the device
timeToDelay = obj.options.calibratorTypeSpecificParamsStruct.LEDWarmupDurationSeconds;
fprintf('Timer will count %2.1f seconds for warming up \n',timeToDelay);
for tt = 1:timeToDelay
    pause(1); % Display a dot per each second to see if the timer is working.
end
disp('Close the timer and the measurement will begin!');

% Then measure it.
%
% The Spectrum range part should be substituted as a function.
measurement = MeasureSpectroradiometer('measurementOption',true);
S = [380 2 201];

end