function setDisplaysInitialState(obj, userPrompt)

    if (obj.options.verbosity > 9)
        fprintf('In PsychImaging.setDisplayInitialState()\n');
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
        obj.nSubprimaries = obj.options.calibratorTypeSpecificParamsStruct.nSubprimaries;
        obj.nPrimaries = obj.options.calibratorTypeSpecificParamsStruct.nPrimaries;
        obj.normalMode = obj.options.calibratorTypeSpecificParamsStruct.normalMode;
        obj.logicalToPhysical = obj.options.calibratorTypeSpecificParamsStruct.logicalToPhysical;
    end
    
    % Specify stereo mode 10 for synchronized flips between left/right displays
    stereoMode = []; % 10; 
    
    % Following for opening a full-screen window
    screenRect = []; 
    
    % Specify pixelSize (30 for 10-bit color, 24 for 8-bit color)
    pixelSize = 24;
    
    % Set background settings.
    backgroundSettings = [1 1 1];    

    %% Display background settings.
    [obj.masterWindowPtr, obj.screenRect] = OpenPlainScreen(backgroundSettings,'screenNum',calStruct.describe.whichScreen-1);
    LoadIdentityClut(obj.masterWindowPtr);
    
    % Blank option for the other display.
    if calStruct.describe.blankOtherScreen
    [obj.slaveWindowPtr, ~] = ...
        OpenPlainScreen(calStruct.describe.blankSettings,'screenNum',calStruct.describe.whichBlankScreen-1);  
         LoadIdentityClut(obj.slaveWindowPtr);
    Screen('Flip', obj.slaveWindowPtr);
    end
   
    % white square for user to focus the spectro-radiometer
    targetSettings = [1 1 1];
    obj.updateBackgroundAndTarget(calStruct.describe.bgColor, targetSettings, calStruct.describe.useBitsPP)    
    
    % Set channel settings here. Just turn all channels on.
    initialChannelSettings = ones(obj.nSubprimaries,obj.nPrimaries); 
    SetChannelSettings(initialChannelSettings);
    
    % Set channel settings here. Just fair channel settings for each screen primary.
%     initialChannelSettings = zeros(obj.nSubprimaries,obj.nPrimaries); 
%     idxChannelTurnOn = [11 12 13 15 16 ; 5 7 8 9 10 ; 1 1 2 3 4]; 
%     for pp = 1:obj.nPrimaries
%         for ii = 1:size(idxChannelTurnOn,2)
%             initialChannelSettings(idxChannelTurnOn(pp,ii),pp) = 1;
%         end
%     end
%     SetChannelSettings(initialChannelSettings);
    
    % Wait for user
    if (userPrompt)
        fprintf('\nHit enter when ready ...');
        FlushEvents;
        GetChar;
        if strcmp(class(obj.radiometerObj), 'SpectroCALdev')
            obj.radiometerObj.switchLaserState(0);
        end
        fprintf('\n\n-------------------------------------------\n');
        fprintf('\nPausing for %d seconds ...', calStruct.describe.leaveRoomTime);
        WaitSecs(calStruct.describe.leaveRoomTime);
        fprintf(' done\n');
        fprintf('\n-------------------------------------------\n\n');
    end
    
end