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
peaks = [422,448,476,474,506,402,532,552,558,592,610,618,632,418,658,632];
index = [1:1:length(peaks)];
whichChannelPrimary1 = 15;
whichChannelPrimary2 = [];
whichChannelPrimary3 = [];

channelIntensityPrimary1 = 1;
channelIntensityPrimary2 = 1;
channelIntensityPrimary3 = 1;

channelSettings = zeros(nChannels, nPrimaries);
channelSettings(whichChannelPrimary1, 1) = channelIntensityPrimary1;
channelSettings(whichChannelPrimary2, 2) = channelIntensityPrimary2;
channelSettings(whichChannelPrimary3, 3) = channelIntensityPrimary3;

% Set channel setting here.
SetChannelSettings(channelSettings);
GetChannelSettings;

%% Which spatial frequency to use.
sineFreqCyclesPerDegOptions = [3, 6, 9, 12, 18];
% while 1
%     inputMessageSpatialFrequency = 'Which spatial frequency to test [3,6,9,12,18]: ';
%     sineFreqCyclesPerDeg = input(inputMessageSpatialFrequency);
%     sineFreqCyclesPerDegOptions = [3, 6, 9, 12, 18];
%
%     if ismember(sineFreqCyclesPerDeg, sineFreqCyclesPerDegOptions)
%         break
%     end
%
%     disp('Spatial frequency should be within the above range!');
% end

% Make a loop here so that we can display all spatial frequencies one by
% one.
nSFs = length(sineFreqCyclesPerDegOptions);
for ss = 1:nSFs
    % Set spatial frequency.
    sineFreqCyclesPerDeg = sineFreqCyclesPerDegOptions(ss);
    
    % Get a key press to start to display image.
    fprintf('Press any key to display image \n',whichChannelPrimary1,sineFreqCyclesPerDeg);
    pause;
    
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
            cyclesPerDeg = sineFreqCyclesPerDeg;
            projectorOption = 'SACCSFA';
            
            if strcmp(projectorOption,'Raw')
                measurementposition = 3;
            end
            
            % On the optical system (SACCSFA, new projector).
            if strcmp(projectorOption,'SACCSFA')
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
                
                % Raw (without SACCSFA, old projector).
            elseif strcmp(projectorOption,'Raw')
                switch cyclesPerDeg
                    case 3
                        barWidthPixel = 80;
                    case 6
                        barWidthPixel = 42;
                    case 9
                        barWidthPixel = 28;
                    case 12
                        barWidthPixel = 21;
                    case 18
                        switch measurementposition
                            case 1
                                barWidthPixel = 12;
                            case 2
                                barWidthPixel = 13;
                                % Nominal
                            case 3
                                barWidthPixel = 14;
                            case 4
                                barWidthPixel = 16;
                            case 5
                                barWidthPixel = 18;
                                % 1.66 d (as of 0814).
                            case 11
                                barWidthPixel = 21;
                                % 2.0 d (as of 0814).
                            case 12
                                barWidthPixel = 26;
                            otherwise
                                barWidthPixel = 14;
                        end
                end
            end
            
            % Calculate the cycles/pixel.
            CalCyclesPerPixel = true;
            if (CalCyclesPerPixel)
                % Resoultion of DLP and camera.
                DLPResolution = [1920 1024];
                CameraResolution = [3088 2064];
                
                % Get the number of cycles on the DLP.
                pixelDLPHorizontal = DLPResolution(1);
                pixelOneCycleOnDLP = barWidthPixel*2;
                numCyclesOnDLP = pixelDLPHorizontal/pixelOneCycleOnDLP;
                
                % Calculate the cycles per pixel on DLP.
                numCyclesOnDLPPerPixel = round(numCyclesOnDLP/pixelDLPHorizontal,3);
            end
            
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
    
    % Print out which image is displaying now.
    fprintf('\t Now displaying - (Ch %d) / (%d cpd) \n',whichChannelPrimary1,sineFreqCyclesPerDeg);
    
end

%% Activate the following commands if you want to quit the screen with key stroke
KEYPRESSTOCLOSE = false;
if (KEYPRESSTOCLOSE)
    KbStrokeWait;
    CloseScreen;
end
