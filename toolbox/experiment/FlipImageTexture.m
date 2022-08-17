function [flipTime] = FlipImageTexture(imageTexture, window, imageWindowRect, options)
% Flip the screen with given PTB texture.
%
% Syntax:
%    [flipTime] = FlipImageTexture(imageTexture, window, imageWindowRect)
%
% Description:
%    This is to display the image by flipping the screen with the given PTB
%    texture.
%
% Inputs:
%    imageTexture               - PTB texture that saves the image.
%    window                     - PTB window for opened screen.
%    imageWindowRect            - Rect corresonding to image texture size.
%
% Outputs:
%    flipTime                   - System flip time of displaying images.
%
% Optional key/value pairs:
%    preFlipTimeDelay           - Default to 0. Make a time delay before
%                                 the flip of the screen. Unit in seconds.
%    afterFlipTimeDelay         - Default to 0. Make a time delay after the
%                                 flip of the screen. Unit in seconds.
%    verbose                    - Boolean. Default true.  Controls plotting
%                                 and printout.
% See also:
%    SetScreenImage, MakeImageTexture

% History:
%    08/17/22      smo          - Wrote it.

%% Set parameters.
arguments
    imageTexture (1,1)
    window (1,1)
    imageWindowRect (1,4)
    options.preFlipTimeDelay (1,1) = 0
    options.afterFlipTimeDelay (1,1) = 0
    options.verbose (1,1) = true
end

%% Flip screen and display the image here.
%
% Draw a texture here.
Screen('DrawTexture', window, imageTexture, [], imageWindowRect);

%% Make pre-flip time delay if you want.
%
% You can make a time delay before and after the flip of the screen. The
% time is corrected using the function not to break the frame.
if (~isempty(options.preFlipTimeDelay))
    preFlipTimeDelay = MatchScreenFrameTime(options.preFlipTimeDelay);
    WaitSecs(preFlipTimeDelay);
end

%% Flip happens here.
flipTime = Screen('Flip', window);

% % Uncomment to test raw speed using flip
% nFrames = 100;
% testTimes = zeros(nFrames,1);
% for ii = 1:nFrames
%     testTimes(ii) = Screen('Flip',window);
% end
% hist(diff(testTimes));

% Show the verbose message if you want.
if (options.verbose)
    fprintf('Image is now being displayed on the screen...\n');
end

%% Make after-flip time delay if you want.
if (~isempty(options.afterFlipTimeDelay))
    afterFlipTimeDelay = MatchScreenFrameTime(options.afterFlipTimeDelay);
    WaitSecs(afterFlipTimeDelay);
end

end
