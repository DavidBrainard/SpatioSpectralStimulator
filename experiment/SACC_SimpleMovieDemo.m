function SimpleMovieDemo(moviename, windowrect)
% Most simplistic demo on how to play a movie.
%
% SimpleMovieDemo(moviename [, windowrect=[]]);
%
% This bare-bones demo plays a single movie whose name has to be provided -
% including the full filesystem path to the movie - exactly once, then
% exits. This is the most minimalistic way of doing it. For a more complex
% demo see PlayMoviesDemo. The remaining demos show more advanced concepts
% like proper timing etc.
%
% The demo will play our standard DualDiscs.mov movie if the 'moviename' is
% omitted.
%

% History:
% 02/05/2009  Created. (MK)
% 06/17/2013  Cleaned up. (MK)

% Check if Psychtoolbox is properly installed:
AssertOpenGL;

if nargin < 1 || isempty(moviename)
    % No moviename given: Use our default movie:
    moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
end

if nargin < 2 || isempty(windowrect)
    windowrect = [];
end

% Wait until user releases keys on keyboard:
KbReleaseWait;

% Select screen for display of movie:
screenid = max(Screen('Screens'));

try
    % Open 'windowrect' sized window on screen, with black [0] background color:
    win = Screen('OpenWindow', screenid, 0, windowrect);
    
    % Open movie file:
    movie = Screen('OpenMovie', win, moviename);
    
    % Start playback engine:
    Screen('PlayMovie', movie, 1);
    
    % Playback loop: Runs until end of movie or keypress:
    index = 1;
    flipTime = zeros(62,1);
    flipTime1 = zeros(62,1);
    texTime = zeros(62,1);
    movieTime = zeros(62,1);
    closeTime = zeros(62,1);
    for ii = 1:62
        tex{ii} = Screen('GetMovieImage', win, movie);
    end

    % Wait until all keys picked up
    while KbCheck; end
    
    for ii = 1:62
        % Wait for next movie frame, retrieve texture handle to it
        startT = GetSecs;
        movieTime(index) = GetSecs-startT;
        
%         % Valid texture returned? A negative value means end of movie reached:
%         if tex<=0
%             % We're done, break out of loop:
%             break;
%         end
        
        % Draw the new texture immediately to screen:
        startT = GetSecs;
        Screen('DrawTexture', win, tex{ii});
        texTime(index) = GetSecs-startT;
        
        % Update display:
        startT = GetSecs;
        flipTime(index) = Screen('Flip', win);
        flipTime1(index) = GetSecs-startT;
        
        % Release texture:
        startT = GetSecs;
        Screen('Close', tex{ii});
        closeTime(index) = GetSecs-startT;
        index = index+1;
        
        % Check the state of the keyboard.
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        
        % If the user is pressing a key, then display its code number and name.
        if keyIsDown
            keyCode = find(keyCode, 1);

            % Note that we use find(keyCode) because keyCode is an array.
            % See 'help KbCheck'
            fprintf('You pressed key %i which is %s\n', keyCode, KbName(keyCode));
            
            % If the user holds down a key, KbCheck will report multiple events.
            % To condense multiple 'keyDown' events into a single event, we wait until all
            % keys have been released.
            KbReleaseWait;
        end
        
    end
    
    figure; hist(diff(flipTime));
    
    % Stop playback:
    Screen('PlayMovie', movie, 0);
    
    % Close movie:
    Screen('CloseMovie', movie);
    
    % Close Screen, we're done:
    sca;
    
catch %#ok<CTCH>
    sca;
    psychrethrow(psychlasterror);
end

return
