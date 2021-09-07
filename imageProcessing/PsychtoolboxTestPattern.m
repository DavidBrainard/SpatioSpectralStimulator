%% Initialize
% tbUse('BrainardLabBase')

%% Clear the workspace and the screen
sca;
close all;
clear all;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer. For help see: Screen Screens?
screens = Screen('Screens');    

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the
% external screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
screenNumber = max(screens);

% Define black and white (white will be 1 and black 0). This is because0
% luminace values are (in general) defined between 0 and 1.onitor to a DisplayPort?

% For help see: help WhiteIndex and help BlackIndex
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window and color it black.
% For help see: Screen OpenWindow?
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

% Get the size of the on screen window in pixels.
% For help see: Screen WindowSize?
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%% Draw lines for Rectangle diagoanl elements
% 
% % 0. So, "full red" is [1 0 0]. "Full green" [0 1 0] and "full blue" [0 0
% % 1]. Play around with these numbers and see the result.
% lineColor = [0 0 0];
% 
% fromH = 0; % starting point of the line (x)
% fromV = 0; % starting point of the line (y)
% toH = screenXpixels; % end point of the line (x)
% toV = screenYpixels; % end point of the line (y)
% 
% fromH2 = screenXpixels; % starting point of the line (x)
% fromV2 = 0; % starting point of the line (y)
% toH2 = 0; % end point of the line (x)
% toV2 = screenYpixels; % end point of the line (y)
% 
% % Screen('DrawDots', window, [dotXpos dotYpos]theClut, dotSizePix, dotColor, [], 2);
% Screen('DrawLine', window, lineColor, fromH, fromV, toH, toV,3); %Screen(‘DrawLine’, windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);
% Screen('DrawLine', window, lineColor, fromH2, fromV2, toH2, toV2,3); %Screen(‘DrawLine’, windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);
% 
% Screen('DrawLine', window, lineColor, 0, 0, screenXpixels, 0,7); % edge horizontal lines
% Screen('DrawLine', window, lineColor, 0, screenYpixels, screenXpixels, screenYpixels,7);
% 
% Screen('DrawLine', window, lineColor, screenXpixels, 0, screenXpixels,screenYpixels,7); % edge vertical lines
% Screen('DrawLine', window, lineColor, 0, 0, 0,screenYpixels,7);

%% Draw cicle of oval shape

% % Take a center of the screen
% [xCenter, yCenter] = RectCenter(windowRect);
% size_circle = 700;
% baseRect = [0 0 size_circle size_circle];
% maxDiameter = max(baseRect) * 1.01;% For Ovals we set a miximum dcodeiameter up to which it is perfect for
% centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter); % Center the rectangle on the centre of the screen
% rectColor = [0 0 1];
% Screen('FillOval', window, rectColor, centeredRect, maxDiameter);
% 
% % Get an displayed image
% % imageArray=Screen('GetImage',window);
% % imwrite(imageArray,'X_image.jpg');
% 
% % Flip to the screen. This command basically draws all of our previous
% % commands onto the screen. For help see: Screen Flip?
% Screen('Flip', window);
% 
% % Now we have drawn to the screen we wait for a keyboard button press (any
% % key) to terminate the demo. For help see: help KbStrokeWaitonitor to a DisplayPort?
% 
% KbStrokeWait;
% 
% % Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
% % all features related to PTB. Note: we leave the variables in the
% % workspace so you can have a look at them if you want.
% % For help see: help sca
% sca;

%% Draw contrast pattern (either vertical or horizontal)

[xCenter, yCenter] = RectCenter(windowRect);
rectColor = [0 0 0];

% Set the barwidth
rectColor = [0 0 0];
barwidth=0; % in pixel unit
logic=2; % 1=vertical / 2=horizontal

% Single rectangle
% baseRect = [0 0 barwidth screenYpixels];
% centeredRect = CenterRectOnPointd(baseRect, xCenter-barwidth/2, yCenter);
% Screen('FillRect', window, rectColor, centeredRect); % Fill in black for the area

% Stripe pattern (black and white)    
if logic==2;
    baseRect = [0 0 barwidth screenYpixels]; % [screenXpixels, screenYpixels] 

    for i=1:1000
        centeredRect_R = CenterRectOnPointd(baseRect, xCenter+barwidth*2*i, yCenter);
        centeredRect_L = CenterRectOnPointd(baseRect, xCenter-barwidth*2*(i-1), yCenter);
        Screen('FillRect', window, rectColor, centeredRect_R);
        Screen('FillRect', window, rectColor, centeredRect_L);
    end

else 
    baseRect = [0 0 screenXpixels barwidth]; % [screenXpixels, screenYpixels] 

    for i=1:1000
        centeredRect_R = CenterRectOnPointd(baseRect, xCenter, yCenter+barwidth*2*i);
        centeredRect_L = CenterRectOnPointd(baseRect, xCenter, yCenter-barwidth*2*(i-1));
        Screen('FillRect', window, rectColor, centeredRect_R);
        Screen('FillRect', window, rectColor, centeredRect_L);
    end
end 

Screen('Flip', window);

%% Activate the following commands if you want to quit the screen with key stroke
% KbStrokeWait;
% sca;