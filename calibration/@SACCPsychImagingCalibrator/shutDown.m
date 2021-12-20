% Method to shutdown the device
function obj = shutDown(obj)
    if (obj.options.verbosity > 9)
        fprintf('In PsychImaging.shutDown() method\n');
    end

    % Close everything.
    CloseSpectroradiometer;
    CloseScreen;
end