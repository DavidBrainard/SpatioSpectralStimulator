% Method to update the stimulus and conduct a single radiometric measurement by
% calling the corresponding method of the attached @Radiometer object.
function [measurement, S] = updateStimulusAndMeasure(obj, bgSettings, targetSettings, useBitsPP)

% Display current subprimary settings.
if (obj.options.verbosity > 1)
    for i=1:obj.nSubprimaries
        fprintf('        Target settings %2.0f   : %2.3f \n\n',i,round((obj.nInputLevels-1)*targetSettings(i)));
    end
end

% Update background and target stimuli.
obj.updateBackgroundAndTarget(bgSettings, targetSettings, useBitsPP);

% Make a delay before the measurement for warming-up the device
timeToDelay = obj.options.calibratorTypeSpecificParamsStruct.LEDWarmupDurationSeconds;
fprintf('        Timer will count %2.1f seconds for warming up \n\n',timeToDelay);
timerForWarmingup = timer('TimerFcn','stat=false','StartDelay',timeToDelay); % Timer setting
start(timerForWarmingup); % Timer starts.
for tt = 1:timeToDelay
    disp('.'); 
    pause(1); % Display a dot per each second to see if the timer is working.
end
delete(timerForWarmingup); % Timer ends.
disp('        Close the Timer and the measurement will begin');

% Then measure it.
measurement = MeasureSPD;
S = WlsToS(measurement);

%     obj.radiometerObj.measure();
%     measurement = obj.radiometerObj.measurement.energy;
%     S = WlsToS(obj.radiometerObj.measurement.spectralAxis(:))
end