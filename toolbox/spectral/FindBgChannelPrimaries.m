function [bgPrimaries,bgSpd,bgCoords] = FindBgChannelPrimaries(targetBgCoords,T,subPrimaryCalstructData, ...
    B_natural,projectIndices,primaryHeadroom,targetLambda,options)
% Find primaries for a desired background spectrum by searching the target background
%
% Syntax:
%    [bgPrimaries,bgSpd,bgCoords] = FindBgChannelPrimaries(targetBgCoords,T,subPrimaryCalstructData, ...
%        B_natural,projectIndices,primaryHeadroom,targetLambda)
%
% Description:
%    This function searches for the desired background SPD which matches in the xy chromaticty
%    based on the initial guess from the input variable 'targetBgxy'.
%    The shape of the final result of the desired SPD is dependent on the variable 'targetLambda'
%    which constraints the smootheness of the curve.
%
% Inputs:
%    targetBgCoords -             The desired color coordinates of the background
%                                 in the color space defined by passed T.
%    T -                          Color matching functions in standard PTB format.
%    subPrimaryCalstructData -    The calibration object that describes the
%                                 device we're working with. It's
%                                 wavelength support needs to match T and
%                                 B_natural.
%    B_natural -                  Basis functions to enforce the smoothness
%                                 on the desired background spectrum.
%                                 This also should be in standard PTB
%                                 format and on the same wavelength support
%                                 as T and the calibration.
%    projectIndices -             Define wavelength range that will be used
%                                 to enforce the smoothnes thorugh the
%                                 projection onto an underlying basis set.
%                                 These are the indices into the wavelength
%                                 vector of the desired range.
%    primaryHeadroom -            Set the headroom to prevent running into
%                                 numerical error at the edges.
%    targetLambda -               This number decides the level of smoothness of
%                                 the final calculation result of the
%                                 desired background spectrum. Higher value
%                                 makes the curve more smoother, but it can
%                                 possibly lose accuracy in matching the
%                                 desired coordinates. You may want to play
%                                 with targetLambda a bit by hand to
%                                 determine tradeoff between smoothness and
%                                 obtaining desired chromaticity.
%
% Outputs:
%    bgPrimaries -                The primaries obtained.
%    bgSpd -                      Spd corresponding to primaries.
%    bgCoords -                   Coordinates corresponding to obtained primaries/spd.
%
% Optional key/value pairs:
%    'ScaleFactor' -              Factor re 1 to determine target
%                                 background luminance, when Scale is true.
%                                 Default 0.5.
%    'Scale' -                    Boolean (default true). Scale target
%                                 coordinates into range near middle of
%                                 projector gamut.
%    'Verbose' -                  Boolean (default true). Plots and
%                                 printouts.
%
% History:
%    09/28/21  smo  Started on it
%    10/01/21  smo, gka, dhb     Clean and vet.

%% Set parameters
arguments
    targetBgCoords
    T
    subPrimaryCalstructData
    B_natural
    projectIndices
    primaryHeadroom
    targetLambda
    options.scaleFactor (1,1) = 0.5
    options.Scale (1,1) = true
    options.Verbose (1,1) = true
end

%% Preparations for the calculations.
%
% This part prepares some matrix and variables which are needed for following calculation parts.
%
% Some transformation matrices.
nCoords = size(T,1);
nBases = size(B_natural,2);
if (nCoords > nBases)
    error('Must pass constraining basis whose dimensionality exceeds that of target color space');
end
M_NaturalToCoords = T*B_natural(:,1:nCoords);
M_CoordsToNatural = inv(M_NaturalToCoords);

% Half-on spectrum.
nSubprimaries = subPrimaryCalstructData.get('nDevices');
halfOnSubprimaries = options.scaleFactor*ones(nSubprimaries,1);
halfOnSpd = PrimaryToSpd(subPrimaryCalstructData,halfOnSubprimaries); % This will be used to set the desired background Y value
halfOnCoords = T*halfOnSpd;

% Find scale of passed coordinates that brings
% them into alignment with the half on coordinates.
if (options.Scale)
    halfOnTargetScaleFactor = targetBgCoords\halfOnCoords;
    targetBgCoordsScale = halfOnTargetScaleFactor*targetBgCoords;
end

% Spectrum range.
S = subPrimaryCalstructData.get('S');
wls = SToWls(S); % This will be used to plot the results after getting the desired background calculation results.

%% Get an initial desired background spd within the appropiate dimension subspace of the natural basis.
desiredBgSpd = B_natural(:,1:nCoords)*M_CoordsToNatural*targetBgCoordsScale;
desiredBgSpd = sum(halfOnSpd)*desiredBgSpd/sum(desiredBgSpd);
desiredBgCoords = T*desiredBgSpd;

%% Search for desired background.
%
% Make a guess at reasonable primary values to start at.
startingBgPrimaries = SpdToPrimary(subPrimaryCalstructData,desiredBgSpd);
startingBgSpd = PrimaryToSpd(subPrimaryCalstructData,startingBgPrimaries);
startingBgCoords = T*startingBgSpd;

% Optimize a desired background spectrum from the above initial guess.
[bgPrimariesIncr] = ReceptorIsolateSpectral(T,desiredBgCoords,subPrimaryCalstructData.get('P_device'),startingBgPrimaries,startingBgPrimaries, ...
    primaryHeadroom,B_natural,projectIndices,targetLambda,subPrimaryCalstructData.get('P_ambient'),'EXCITATIONS',true);
bgPrimaries = startingBgPrimaries + bgPrimariesIncr;
bgSettings = PrimaryToSettings(subPrimaryCalstructData,bgPrimaries);
bgPrimariesQuantized = SettingsToPrimary(subPrimaryCalstructData,bgSettings);
bgSpd = PrimaryToSpd(subPrimaryCalstructData,bgPrimaries);
bgCoords = T*bgSpd;

% Plot the results of the starting BG Spd and the obtained Bg spd.
if (options.Verbose)
    figure; clf; hold on;
    plot(wls,startingBgSpd,'r','LineWidth',4);
    plot(wls,bgSpd,'k','LineWidth',2);
    legend('Starting BG','Obtained BG');
    xlabel('Wavelength (nm)');
    ylabel('Spectral power');

    % Print out mean value of background primaries
    fprintf('        Mean value of background primaries: %0.2f\n',mean(bgPrimaries));

    % Print out how we did at getting the desired background coordinates
    for cc = 1:nCoords
        fprintf('        Desired   coorindate %d: %0.3g\n',cc,desiredBgCoords(cc));
        fprintf('        Starting  coorindate %d: %0.3g\n',cc,startingBgCoords(cc));
        fprintf('        Obtained  coorindate %d: %0.3g\n',cc,bgCoords(cc));
    end
end

end
