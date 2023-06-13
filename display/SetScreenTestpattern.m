% SetScreenTestpattern
% 
% It displays some test patterns as desired on the screen. It is able to
% display contrast pattern, circle on the plain screen, etc.

% History:
%    12/06/21 smo    - Cleaned up.
%    09/14/22 smo    - Added an option to control projector LED channel
%                      settings.

%% Start over.
sca;
close all;
clear all;

%% Open the projector.
initialScreenSettings = [1 1 1];
[window, windowRect] = OpenPlainScreen(initialScreenSettings);
screenXpixels = windowRect(3);
screenYpixels = windowRect(4);

%% Set the projector LED channel settings.
nPrimaries = 3;
nChannels = 16;

% Set which channel to use per each primary. We will use only two
% primaries.
%
% Peak wavelength in order of channel number. Note that it is not ascending
% order.
% [422,448,476,474,506,402,532,552,558,592,610,618,632,418,658,632]
peaks = [422,448,476,474,506,402,532,552,558,592,610,618,632,418,658,632];
index = [1:1:length(peaks)];
whichChannelPrimary1 = 15;
whichChannelPrimary2 = [];
whichChannelPrimary3 = [];

channelIntensityPrimary1 = 1;
channelIntensityPrimary2 = 0.3;
channelIntensityPrimary3 = 0.3;

channelSettings = zeros(nChannels, nPrimaries);
channelSettings(whichChannelPrimary1, 1) = channelIntensityPrimary1;
channelSettings(whichChannelPrimary2, 2) = channelIntensityPrimary2;
channelSettings(whichChannelPrimary3, 3) = channelIntensityPrimary3;

% Set channel setting here.
SetChannelSettings(channelSettings);
GetChannelSettings;

%% Set which type of screen to display here.
% 
% Set one among [Circle, Contrast, SingleLine]
ScreenPatternType = 'Contrast';

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
        % 18 cpd = 14 pixels / 12 cpd = 21 pixels / 9 cpd = 28 pixels / 6
        % cpd = 42 pixels / 3 cpd = 80 pixels.
        rectColor = [0 0 0];
        
        cyclesPerDeg = 3;
        
        % On the optical system.
        switch cyclesPerDeg
            case 3
                barWidthPixel = 23;
            case 6
                barWidthPixel = 11;
            case 9
                barWidthPixel = 8;
            case 12
                barWidthPixel = 6;
            case 18
                barWidthPixel = 4;
        end
        
        % Raw (old projector).
%         switch cyclesPerDeg
%             case 3
%                 barWidthPixel = 80;
%             case 6
%                 barWidthPixel = 42;
%             case 9
%                 barWidthPixel = 28;
%             case 12
%                 barWidthPixel = 21;
%             case 18
%                 barWidthPixel = 14;
%         end
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
KEYPRESSTOCLOSE = false;
if (KEYPRESSTOCLOSE)
    KbStrokeWait;
    CloseScreen;
end
