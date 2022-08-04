function [flipTime] = DisplayScreenPattern(window,windowRect,options)
% Display a pattern on the screen using PTB.
%
% Syntax:
%    [] = DisplayScreenPattern(window,windowRect)
%
% Description:
%    This displays a pattern using PTB. It can display a circle, a single
%    line, a contrast pattern, and a cross bar at the center of the screen.
%    This will be useful to add a slide between
%
% Inputs:
%    window -                     PTB window for opened screen.
%    windowRect -                 Rect corresonding to window.
%
% Outputs:
%    flipTime -                   System time of screen flip in seconds.
%
% Optional key/value pairs:
%    patternType -                Default to 'crossbar'. Choose a desired
%                                 pattern type to display within four.
%                                 [circle, contrast, singleLine, crossbar]
%    patternColor -               Deafult to white ([1 1 1]). Set color of
%                                 the pattern. Each value should be ranged
%                                 within 0-1.
%    imageBackground -            Image to use as a background instead of
%                                 plain screen. Input should be in format
%                                 either double or uint8.
%    timeDelay -                  Default to 0. Make a time delay before
%                                 making a flip of the image.
%    verbose -                    Boolean. Default true. Controls
%                                 printout.
%
% See also:
%    N/A

% History:
%   02/10/22 smo                  Started on it.
%   02/15/22 smo                  Now it is possible to display crossbar
%                                 pattern on image.
%   08/03/22 smo                  Added an option to make a time delay
%                                 before displaying an image.

%% Set parameters.
arguments
    window (1,1)
    windowRect (1,4)
    options.patternColor (1,3) = [1 1 1]
    options.patternType = 'crossbar'
    options.imageBackground = []
    options.preFlipTimeDelay (1,1) = 0
    options.afterFlipTimeDelay (1,1) = 0
    options.verbose (1,1) = true
end

%% Get screen info.
screenXpixels = windowRect(3);
screenYpixels = windowRect(4);

%% Display a pattern here.
switch options.patternType
    case 'circle'
        % We will place a rectangle at the center of the screen.
        [xCenter, yCenter] = RectCenter(windowRect);
        sizeCirclePixel = 200;
        baseRect = [0 0 sizeCirclePixel sizeCirclePixel];
        maxDiameter = max(baseRect) * 1.01;
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
        Screen('FillOval', window, options.patternColor, centeredRect, maxDiameter);
        Screen('Flip', window);
        
    case 'contrast'
        [xCenter, yCenter] = RectCenter(windowRect);
        
        % Set the barwidth and directions.
        barWidthPixel = 20;
        whichSideBar = 'vertical';
        
        % Set stripe pattern.
        switch whichSideBar
            case 'horizontal'
                baseRect = [0 0 screenXpixels barWidthPixel];
                for i=1:1000
                    centeredRect_R = CenterRectOnPointd(baseRect, xCenter, yCenter+barWidthPixel*2*i);
                    centeredRect_L = CenterRectOnPointd(baseRect, xCenter, yCenter-barWidthPixel*2*(i-1));
                    Screen('FillRect', window, options.patternColor, centeredRect_R);
                    Screen('FillRect', window, options.patternColor, centeredRect_L);
                end
                
            case 'vertical'
                baseRect = [0 0 barWidthPixel screenYpixels];
                for i=1:1000
                    centeredRect_R = CenterRectOnPointd(baseRect, xCenter+barWidthPixel*2*i, yCenter);
                    centeredRect_L = CenterRectOnPointd(baseRect, xCenter-barWidthPixel*2*(i-1), yCenter);
                    Screen('FillRect', window, options.patternColor, centeredRect_R);
                    Screen('FillRect', window, options.patternColor, centeredRect_L);
                end
        end
        
        Screen('Flip', window);
        
    case 'singleline'
        barWidthPixel = 1;
        whichSideBar = 'horizontal';
        
        switch whichSideBar
            case 'horizontal'
                fromX = 0;
                fromY = 0.5 * screenYpixels;
                toX = screenXpixels;
                toY = 0.5 * screenYpixels;
                Screen('DrawLine',window,options.patternColor,fromX,fromY,toX,toY,barWidthPixel);
                Screen('Flip', window);
                
            case 'vertical'
                fromX = 0.5 * screenXpixels;
                fromY = 0;
                toX = 0.5 * screenXpixels;
                toY = screenYpixels;
                Screen('DrawLine',window,options.patternColor,fromX,fromY,toX,toY,barWidthPixel);
                Screen('Flip', window);
        end
        
    case 'crossbar'
        % Set up alpha-blending for smooth (anti-aliased) lines.
        Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
        
        % Get the centre coordinate of the window.
        [xCenter, yCenter] = RectCenter(windowRect);
        
        % Here we set the size of the arms of our fixation cross
        fixCrossDimPix = 30;
        
        % Now we set the coordinates (these are all relative to zero we
        % will let the drawing routine center the cross in the center of
        % our monitor for us).
        xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
        yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
        allCoords = [xCoords; yCoords];
        
        % Set the line width for our fixation cross.
        lineWidthPix = 2;
        
        % Draw the fixation cross here at the center of the screen.
        nInputLevels = 256;
        patternColor = SettingsToIntegers(options.patternColor,'nInputLevels',nInputLevels);
        
        % Display image as a background if you want.
        if (~isempty(options.imageBackground))
            % We will blend crossbar on the image background.
            % This version assumes the image is in a square shape.
            crossbarImage = options.imageBackground;
            imageSize = size(crossbarImage);
            imageCenter = imageSize(1) * 0.5;
            imageLineWidth = [imageCenter-lineWidthPix : imageCenter+lineWidthPix];
            crossbarCoords = [imageCenter-fixCrossDimPix/2 : imageCenter+fixCrossDimPix/2];
            
            % Set crossbar image here and display.
            for ii = 1:size(options.patternColor,2)
                % Horizontal part of the crossbar.
                crossbarImage(imageLineWidth,crossbarCoords,ii) = options.patternColor(ii);
                % Vertical part of the crossbar.
                crossbarImage(crossbarCoords,imageLineWidth,ii) = options.patternColor(ii);
            end
            [flipTime] = SetScreenImage(crossbarImage, window, windowRect, ...
                'preFlipTimeDelay', options.preFlipTimeDelay, 'afterFlipTimeDelay', options.afterFlipTimeDelay, ...
                'verbose', options.verbose);
            
        else
            % Display crossbar on the plain screen.
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, [xCenter yCenter], 2);
            
            % Flip to the screen.
            Screen('Flip', window);
        end
        
    otherwise
end

if (options.verbose)
    fprintf('(%s) pattern is displaying...\n',options.patternType);
end

end