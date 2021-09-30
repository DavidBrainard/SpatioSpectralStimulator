function [obtainedBgSpd,obtainedBgxyY] = FindDesiredBackgroundPrimaries(T_xyz,T_cones,B_natural, ... 
                                       targetBgxy,subPrimaryCalstructData,subPrimarySpd,ambientSpd,projectIndices,primaryHeadRoom,targetLambda)
% Find a desired background spectrum by searching the target background 
% xy chromaticity values.
% 
% Syntax:
%    [obtainedBgSpd,obtainedBgxyY] = FindDesiredBackgroundPrimaries(T_xyz,T_cones,B_natural,...,
%                                    targetBgxy,subPrimaryCalstructData,subPrimarySpd,ambientSpd,projectIndices,primaryHeadRoom,targetLambda)
% Description:
%    This function searches for the desired background SPD which matches in the xy chromaticty 
%    based on the initial guess from the input variable 'targetBgxy'. 
%    The shape of the final result of the desired SPD is dependent on the variable 'targetLambda' 
%    which constraints the smootheness of the curve. 
%
% Inputs:
%    T_xyz -                      Color matching functions in standard PTB format.
%    T_cones -                    Spectral cone sensitivities in standard PTB format.
%    B_natural -                  Basis functions to enforce the smoothness
%                                 on the desired background spectrum. 
%                                 This also should be in standard PTB
%                                 format. 
%    targetBgxy -                 The CIE xy chromaticities of the target background. 
%                                 This will be used as an initial guess for searching the
%                                 desired background xyY.
%    subPrimaryCalstructData -    The subprimary calstructdata to calculate
%                                 some routines including
%                                 'SpdToPrimary','PrimaryToSPd','SettingsToPrimary',
%                                 and 'PrimaryToSettings'
%                                 Also, the information on the variables,
%                                 'nSubprimaries' and 'wls', in the script, 
%                                 are acquired from here.
%    subPrimarySpd -              Spectral power distribution of all subprimaries.
%    ambientSpd -                 Spectral power distribution of the ambient light.
%    projectIndices -             Define wavelength range that will be used to enforce the smoothnes
%                                 thorugh the projection onto an underlying basis set.
%    PrimaryHeadRoom -            Set the headroom to prevent running into
%                                 numerical error at the edges.
%    targetLambda -               This number decides the level of smoothness of
%                                 the final calculation result of the
%                                 desired background spectrum. Higher value
%                                 makes the curve more smoother, 
%                                 but it can possibly lose accuracy
%                                 in matching the desired xy chromaticity.
%                                 You may want to play with targetLambda a bit 
%                                 by hand to determine tradeoff between
%                                 smoothness and obtaining desired chromaticity.
%
% Outputs:
%    obtainedBgSpd -              Final optimized result of the desired
%                                 background SPD.
%    obtainedBgxyY -              Final optimized result of the desired 
%                                 background xy chromaticity.
%
% Optional key/value pairs:
%    N/A
%
% History: 
%    09/28/21  smo  Started on it

%% Preparations for the calculations.
%
% This part prepares some matrix and variables which are needed for following calculation parts.
%
% Some transformation matrices.
M_NaturalToXYZ = T_xyz*B_natural(:,1:3);
M_XYZToNatural = inv(M_NaturalToXYZ);
M_NaturalToLMS = T_cones*B_natural(:,1:3);
M_LMSToNatural = inv(M_NaturalToLMS);

% Half-on spectrum.
nSubprimaries = subPrimaryCalstructData.get('nDevices');
halfOnSubprimaries = 0.5*ones(nSubprimaries,1);
halfOnSpd = PrimaryToSpd(subPrimaryCalstructData,halfOnSubprimaries); % This will be used to set the desired background Y value

% Spectrum range.
S = subPrimaryCalstructData.get('S');
wls = SToWls(S); % This will be used to plot the results after getting the desired background calculation results.

%% Get the desired background chromaticity.
%
% Start by defining a smooth spectrum with the chromaticity
% we want.  This is not necessarily within the span of the
% subprimaries we have, but serves as a starting point.
startingTargetBgY = 1; % Set an arbitrary number for Y to start.
targetBgxyY = [targetBgxy ; startingTargetBgY];
targetBgXYZ = xyYToXYZ(targetBgxyY);
desiredBgSpd = B_natural(:,1:3)*M_XYZToNatural*targetBgXYZ;
desiredBgSpd = sum(halfOnSpd)*desiredBgSpd/sum(desiredBgSpd);
desiredBgXYZ = T_xyz*desiredBgSpd;
desiredBgxyY = XYZToxyY(desiredBgXYZ); % This gives the same result in xy chromaticity as 'targetBgxyY', but different Y value which was set to sum of half-on spectrums.

%% Search for desired background.
%
% Make a guess at reasonable primary values to start at.
startingBgPrimaries = SpdToPrimary(subPrimaryCalstructData,desiredBgSpd);
startingBgSpd = PrimaryToSpd(subPrimaryCalstructData,startingBgPrimaries);
startingBgXYZ = T_xyz*startingBgSpd;
startingBgxyY = XYZToxyY(startingBgXYZ);

% Optimize a desired background spectrum from the above initial guess.
[bgPrimariesIncr] = ReceptorIsolateSpectral(T_xyz,desiredBgXYZ,subPrimarySpd,startingBgPrimaries,startingBgPrimaries, ...
                    primaryHeadRoom,B_natural,projectIndices,targetLambda,ambientSpd,'EXCITATIONS',true);
bgPrimaries = startingBgPrimaries + bgPrimariesIncr;
bgSettings = PrimaryToSettings(subPrimaryCalstructData,bgPrimaries);
bgPrimariesQuantized = SettingsToPrimary(subPrimaryCalstructData,bgSettings);
obtainedBgSpd = PrimaryToSpd(subPrimaryCalstructData,bgPrimaries);
obtainedBgLMS = T_cones*obtainedBgSpd;
obtainedBgXYZ = T_xyz*obtainedBgSpd;
obtainedBgxyY = XYZToxyY(obtainedBgXYZ);

% Plot the results of the starting BG Spd and the obtained Bg Spd.
figure; clf; hold on;
plot(wls,startingBgSpd,'r','LineWidth',4);
plot(wls,obtainedBgSpd,'k','LineWidth',2);
legend('Starting BG','Obtained BG');
xlabel('Wavelength (nm)');
ylabel('Spectral power');

% Print out the results in the CIE xy chromaticity.
fprintf('        Desired  background x,y = %0.3f,%0.3f\n',desiredBgxyY(1),desiredBgxyY(2));
fprintf('        Starting background x,y = %0.3f,%0.3f\n',startingBgxyY(1),startingBgxyY(2));
fprintf('        Obtained background x,y = %0.3f,%0.3f\n',obtainedBgxyY(1),obtainedBgxyY(2));
fprintf('        Mean value of background primaries: %0.2f\n',mean(bgPrimaries));

end
