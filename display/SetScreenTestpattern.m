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

%% Initialize.
%
% Open the screen.
initialScreenSettings = [1 1 1];
[window, windowRect] = OpenPlainScreen(initialScreenSettings);
screenXpixels = windowRect(3);
screenYpixels = windowRect(4);

% Set the channel settings.
nPrimaries = 3;
nChannels = 16;

channelSettings = ones(nChannels, nPrimaries);
channelIntensity = 0.1;
channelSettings = channelIntensity * channelSettings;

SetChannelSettings(channelSettings);

%% Set which type of screen to display here.
% 
% Set one among [Circle, Contrast, SingleLine]
ScreenPatternType = 'SingleLine';

switch ScreenPatternType
     case 'Circle'
        [xCenter, yCenter] = RectCenter(windowRect);
        sizeCirclePixel = 200;
        baseRect = [0 0 sizeCirclePixel sizeCirclePixel];
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
            case 'horizontal'
            baseRect = [0 0 screenXpixels barWidthPixel]; 
            for i=1:1000
                centeredRect_R = CenterRectOnPointd(baseRect, xCenter, yCenter+barWidthPixel*2*i);
                centeredRect_L = CenterRectOnPointd(baseRect, xCenter, yCenter-barWidthPixel*2*(i-1));
                Screen('FillRect', window, rectColor, centeredRect_R);
                Screen('FillRect', window, rectColor, centeredRect_L);
            end
            
            case 'vertical'
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
        % Set the barwidth and directions.
        lineColor = [0 0 0];
        barWidthPixel = 1; 
        whichSideBar = 'horizontal'; 
        
        switch whichSideBar
            case 'horizontal'
            fromX = 0;
            fromY = 0.5 * screenYpixels;
            toX = screenXpixels;
            toY = 0.5 * screenYpixels;
            Screen('DrawLine',window,lineColor,fromX,fromY,toX,toY,barWidthPixel);
            Screen('Flip', window);

            case 'vertical'
            fromX = 0.5 * screenXpixels;
            fromY = 0;
            toX = 0.5 * screenXpixels;
            toY = screenYpixels;
            Screen('DrawLine',window,lineColor,fromX,fromY,toX,toY,barWidthPixel);       
            Screen('Flip', window);
        end
        
    otherwise 
end

%% Activate the following commands if you want to quit the screen with key stroke
KbStrokeWait;
CloseScreen;
