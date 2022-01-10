function updateBackgroundAndTarget(obj, bgSettings, targetSettings, useBitsPP)
% If the user sets useBitsPP to true, something has gone wrong.
if (useBitsPP)
    error('The PsychImaging calibrator does not support bits++ yet');
end

% Set and check that the number of subprimaries is right
if (length(targetSettings) ~= obj.nSubprimaries)
    error(sprintf('Wrong number %d of targetSettings entries, should be %d',length(targetSettings),obj.nSubprimaries));
end

% Subprimary values should be between 0 and 1.  We evenually
% multiply by 252 round to integer values.
if (any(targetSettings < 0 | targetSettings > 1))
    error('Entries of targetSettings should be between 0 and 252');
end

try
    % Ingore the bgSettings, not meaningful for the channel
    % calibration.
    %
    % Set the channel of whichever primary we're using to the
    % values in targetSettings.
    %
    % Set other two screen primaries to the arbitraryBlack setting, to get
    % ambient measurement out of the mud.
    
    % Check the target screen primary is within the working range.
    if (obj.whichPrimary > obj.nPrimaries)
        error('SACC display has only three screen primaries');
    end
    
    % Set the target subprimary settings here.
    allScreenPrimaries = [1:1:obj.nPrimaries];
    otherScreenPrimaries = setdiff(allScreenPrimaries,obj.whichPrimary);
    
    channelSettings = zeros(obj.nSubprimaries,obj.nPrimaries); % Base matrix for the subprimary settings.
    channelSettings(:,obj.whichPrimary) = targetSettings'; % Target primary setting.
    channelSettings(:,otherScreenPrimaries) = obj.arbitraryBlack; % Other primaries settings. 
    
    SetChannelSettings(channelSettings);
    
catch err
    sca;
    rethrow(err);
end
end
