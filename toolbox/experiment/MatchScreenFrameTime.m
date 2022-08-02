function [correctTime] = MatchScreenFrameTime(time,options)
% Correct the time that that does not break the screen frame rate.
%
% Syntax:
%    [correctTime] = MatchScreenFrameTime(time)
%
% Description:
%    This function returns a corrected time according to display's frame
%    rate. As we don't want to crash the frames between short-term display
%    of the images, so this would make sure that the time is a mulitiple of
%    the 1/frame rate.
%
% Inputs:
%     time                          - Time that you want to correct in
%                                     seconds unit.
%
% Outputs:
%     correctedTime                 - Corrected time that is multiple of 
%                                     1/(frame rate).
%
% Optional key/value pairs:
%    frameRate                      - Default to 119.9573 hz. We use 120 hz
%                                     display for SACC project, whose
%                                     1/(frame rate) is about 0.0083 seconds.
%
% See also
%    N/A

% History:
%    07/22/22  smo                  - Wrote it.

%% Set parameters.
arguments
    time (1,1) 
    options.frameRate (1,1) = 119.9573
end

%% Calculate the corrected time.
ifi = 1/options.frameRate;
nUnitsFrameTime = round(time/ifi);

% Print out the corrected time.
correctTime = ifi * nUnitsFrameTime;

end 