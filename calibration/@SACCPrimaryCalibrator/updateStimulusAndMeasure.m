% Method to update the stimulus and conduct a single radiometric measurement by
% calling the corresponding method of the attached @Radiometer object.
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
fprintf('        Timer will count %2.1f seconds for warming up \n\n',timeToDelay);
for tt = 1:timeToDelay
    disp('.'); 
    pause(1); % Display a dot per each second to see if the timer is working.
end
disp('        Close the timer and the measurement will begin!');

% Then measure it. (THIS PART SHOULD BE MODIFIED LATER)
measurement = MeasureSpectroradiometer('measurementOption',true);
% spectralAxis = linspace(380,780,size(measurement,1))'; % cf. spectralAxis = [380,382,...,780]
% S = WlsToS(spectralAxis);
S = [380 2 201];

% (ORIGINAL CODE)
%     obj.radiometerObj.measure();
%     measurement = obj.radiometerObj.measurement.energy;
%     S = WlsToS(obj.radiometerObj.measurement.spectralAxis(:))
end