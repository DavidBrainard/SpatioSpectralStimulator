function [flipTime] = FlipImageTexture(imageTexture, window, imageWindowRect, options)
% Make a PTB texture of an image.
%
% Syntax:
%    [imageTexture] = MakeImageTexture(image,window,windowRect)
%
% Description:
%    This is to display test images in the experiment for the SACC project
%    using PTB.
%
% Inputs:
%    image -                      Test images to display on the screen.
%                                 This should be in a image format, not a
%                                 cal format. For example, 512 x 512 x 3
%                                 in double.
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Outputs:
%    flipTime                     System flip time of displaying images.
%
% Optional key/value pairs:
%    timeDelay                    Default to 0. Make a time delay before
%                                 making a flip of the image.
%    addFixationPointImage        Default to false. If it is set to true,
%                                 add a fixation point at the center of the
%                                 image. This is useful when you want to
%                                 add one on the stimuli for SACC project.
%    verbose                      Boolean. Default true.  Controls plotting
%                                 and printout.
% See also:
%    SetScreenImage, MakeImageTexture

% History:
%    08/17/22      smo          - Wrote it.

%% Set parameters.
arguments
    imageTexture
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

% You can make a time delay before and after the flip of the screen. The
% time is corrected using the function not to break the frame.
%
% Make pre-flip time delay.
if (~isempty(options.preFlipTimeDelay))
    preFlipTimeDelay = MatchScreenFrameTime(options.preFlipTimeDelay);
    WaitSecs(preFlipTimeDelay);
end

% Flip happens here.
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

% Make after-flip time delay.
if (~isempty(options.afterFlipTimeDelay))
    afterFlipTimeDelay = MatchScreenFrameTime(options.afterFlipTimeDelay);
    WaitSecs(afterFlipTimeDelay);
end

end