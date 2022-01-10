function setDisplaysInitialState(obj, userPrompt)

if (obj.options.verbosity > 9)
    fprintf('In SACCPrimary.setDisplayInitialState()\n');
end

% Make a local copy of obj.cal so we do not keep calling it and regenerating it
calStruct = obj.cal;

%  Get identity clut if useBitsPP is enabled
if ((calStruct.describe.useBitsPP) && (isempty(obj.identityGammaForBitsPP)))
    error('Support for bitsPP has not been implemented in @PsychImagingCalibrator !');
end % if (calStruct.config.useBitsPP)

% Retrieve passed custom params
if (~isempty(obj.options.calibratorTypeSpecificParamsStruct))
    % SACC Subprimary calibration settings
    obj.whichPrimary = obj.options.calibratorTypeSpecificParamsStruct.whichPrimary;   
    obj.nSubprimaries = obj.options.calibratorTypeSpecificParamsStruct.nSubprimaries;
    obj.nInputLevels = obj.options.calibratorTypeSpecificParamsStruct.nInputLevels;
    obj.normalMode = obj.options.calibratorTypeSpecificParamsStruct.normalMode;
    obj.arbitraryBlack = obj.options.calibratorTypeSpecificParamsStruct.arbitraryBlack;
    obj.logicalToPhysical = obj.options.calibratorTypeSpecificParamsStruct.logicalToPhysical;
end

% Disable syncing, we do not care for this kind of calibration (regular
% single screen, not stereo, not Samsung)
Screen('Preference', 'SkipSyncTests', 1);

% Start PsychImaging.
% stereoMode = [];
% screenRect = [];
% pixelSize = 24;
% PsychImaging('PrepareConfiguration');

% Set background settings.
backgroundSettings = [1 1 1];    

%% WORKING ON THIS PART
% Display background settings. (NEW)
[obj.masterWindowPtr, obj.screenRect] = OpenPlainScreen(backgroundSettings,'screenNum',calStruct.describe.whichScreen-1);
LoadIdentityClut(obj.masterWindowPtr);

% (ORIGINAL CODE)
% [obj.masterWindowPtr, obj.screenRect] = ...
%     PsychImaging('OpenWindow', calStruct.describe.whichScreen-1, 255*backgroundSettings, screenRect, pixelSize, [], stereoMode);
% LoadIdentityClut(obj.masterWindowPtr);

% Blank option for the other display. (NEW)
if calStruct.describe.blankOtherScreen
    [obj.slaveWindowPtr, ~] = ...
        OpenPlainScreen(calStruct.describe.blankSettings,'screenNum',calStruct.describe.whichBlankScreen-1);  
         LoadIdentityClut(obj.slaveWindowPtr);
    Screen('Flip', obj.slaveWindowPtr);
end

%  (ORIGINAL CODE)
% if calStruct.describe.blankOtherScreen
%     
%     [obj.slaveWindowPtr, ~] = ...
%         PsychImaging('OpenWindow', calStruct.describe.whichBlankScreen-1, 255*calStruct.describe.blankSettings, [], pixelSize, [], stereoMode);
%     
%     LoadIdentityClut(obj.slaveWindowPtr);
%     Screen('Flip', obj.slaveWindowPtr);
% end

%%
% White square for user to focus the spectro-radiometer
targetSettings = ones(1,obj.nSubprimaries);
obj.updateBackgroundAndTarget(backgroundSettings, targetSettings, calStruct.describe.useBitsPP)

% Check if the number of subprimaris matches.
if (length(obj.logicalToPhysical) ~= obj.nSubprimaries)
    error('Mismatch in number of subprimaries specificaiton and logical to physical array');
end

% Check the number of primary to calibrate.
if (obj.whichPrimary > obj.nPrimaries) 
    error('SACC display has only three primaries');
end

% Set subprimary settings here. Projector target primary is turn at all
% subprimary settings, while the other primaries are turned off.
initialChannelSettings = zeros(obj.nSubprimaries,obj.nPrimaries); % Base matrix for subprimary settings.
initialChannelSettings(:,obj.whichPrimary) = 1;
SetChannelSettings(initialChannelSettings);

% Wait for user
if (userPrompt)
    fprintf('\nHit enter when ready ...');
    FlushEvents;
    GetChar;
    fprintf('\n\n-------------------------------------------\n');
    fprintf('\nPausing for %d seconds ...', calStruct.describe.leaveRoomTime);
    WaitSecs(calStruct.describe.leaveRoomTime);
    fprintf(' done\n');
    fprintf('\n-------------------------------------------\n\n');
end
end