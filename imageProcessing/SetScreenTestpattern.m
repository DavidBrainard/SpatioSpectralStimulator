% SetScreenTestpattern
% 
% It displays some test patterns as desired on the screen. It is able to
% display contrast pattern, circle on the plain screen, etc.

% History:
%    12/6/21 smo   Cleaned up.
%% Start over.
sca;
close all;
clear all;

%% Set basic Psychtoolbox settings.
%
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
screens = Screen('Screens');    
screenNumber = max(screens);

% Define black and white (white will be 1 and black 0). 
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open window in either white or black.
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

% Get the size of the on screen window in pixels.
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Set which type of screen to display here.
% 
% Set one among [Rectangle, Circle, Contrast, SingleLine]
ScreenPatternType = 'SingleLine';

switch ScreenPatternType
    case 'Rectangle'
        
        % Draw lines for Rectangle diagoanl elements
        lineColor = [0 0 0];

        fromH = 0; % starting point of the line (x)
        fromV = 0; % starting point of the line (y)
        toH = screenXpixels; % end point of the line (x)
        toV = screenYpixels; % end point of the line (y)

        fromH2 = screenXpixels; % starting point of the line (x)
        fromV2 = 0; % starting point of the line (y)
        toH2 = 0; % end point of the line (x)
        toV2 = screenYpixels; % end point of the line (y)

        Screen('DrawDots', window, [dotXpos dotYpos], theClut, dotSizePix, dotColor, [], 2);
        Screen('DrawLine', window, lineColor, fromH, fromV, toH, toV,3); %Screen(‘DrawLine’, windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);
        Screen('DrawLine', window, lineColor, fromH2, fromV2, toH2, toV2,3); %Screen(‘DrawLine’, windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);

        Screen('DrawLine', window, lineColor, 0, 0, screenXpixels, 0,7); % edge horizontal lines
        Screen('DrawLine', window, lineColor, 0, screenYpixels, screenXpixels, screenYpixels,7);

        Screen('DrawLine', window, lineColor, screenXpixels, 0, screenXpixels,screenYpixels,7); % edge vertical lines
        Screen('DrawLine', window, lineColor, 0, 0, 0,screenYpixels,7);

    case 'Circle'
        [xCenter, yCenter] = RectCenter(windowRect);
        size_circle = 200;
        baseRect = [0 0 size_circle size_circle];
        maxDiameter = max(baseRect) * 1.01;% For Ovals we set a miximum dcodeiameter up to which it is perfect for
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter); % Center the rectangle on the centre of the screen
        rectColor = [0 0 0];
        Screen('FillOval', window, rectColor, centeredRect, maxDiameter);

        Screen('Flip', window);

    case 'Contrast'
        [xCenter, yCenter] = RectCenter(windowRect);
        rectColor = [0 0 0];

        % Set the barwidth and directions.
        rectColor = [0 0 0];
        barWidthPixel = 20; 
        whichSideBar = 'vertical'; 

        % Set stripe pattern.    
        switch whichSideBar
            case 'vertical'
            baseRect = [0 0 screenXpixels barWidthPixel]; 
            for i=1:1000
                centeredRect_R = CenterRectOnPointd(baseRect, xCenter, yCenter+barWidthPixel*2*i);
                centeredRect_L = CenterRectOnPointd(baseRect, xCenter, yCenter-barWidthPixel*2*(i-1));
                Screen('FillRect', window, rectColor, centeredRect_R);
                Screen('FillRect', window, rectColor, centeredRect_L);
            end
            
        % Horizontal pattern.
            case 'horizontal'
            baseRect = [0 0 barWidthPixel screenYpixels];       
            for i=1:1000
                centeredRect_R = CenterRectOnPointd(baseRect, xCenter+barWidthPixel*2*i, yCenter);
                centeredRect_L = CenterRectOnPointd(baseRect, xCenter-barWidthPixel*2*(i-1), yCenter);
                Screen('FillRect', window, rectColor, centeredRect_R);
                Screen('FillRect', window, rectColor, centeredRect_L);
            end         
        end

        Screen('Flip', window);

    case 'SingleLine'
        [xCenter, yCenter] = RectCenter(windowRect);
        rectColor = [0 0 0];

        % Set the barwidth
        rectColor = [0 0 0];
        barWidthPixel = 5; % in pixel unit
        whichSideBar = 'vertical'; 
        
        % Single rectangle
    switch whichSideBar
        case 'vertical'
        baseRect = [0 0 barWidthPixel screenYpixels];
        centeredRect = CenterRectOnPointd(baseRect, xCenter-barWidthPixel/2, yCenter);
        Screen('FillRect', window, rectColor, centeredRect); % Fill in black for the area
        Screen('Flip', window);
        
        case 'horizontal'
        baseRect = [0 0 screenXpixels barWidthPixel];
        centeredRect = CenterRectOnPointd(baseRect, xCenter-barWidthPixel/2, yCenter);
        Screen('FillRect', window, rectColor, centeredRect); % Fill in black for the area
        Screen('Flip', window);
    end
    otherwise 
end

%% Activate the following commands if you want to quit the screen with key stroke
KbStrokeWait;
sca;
