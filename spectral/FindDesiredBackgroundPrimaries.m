function [bgSpd,bgxyY] = FindDesiredBackgroundPrimaries(targetBgxy,subPrimaryCalData,primaryHeadRoom,targetLambda)
% Find the desired background spectrum starting from the initial background
% setting.
% 
% Syntax:
%   
%
% Description:
%    
%
% Inputs:
%    targetBgxy
%    subPrimaryCalData
%    PrimaryHeadRoom
%    targetLambda
%
% Outputs:
%    bgSpd
%    bgxyY
%
% Optional key/value pairs:
%    
%
% History: 
%    09/28/21  dhb,smo  Started on it

%% Some transformation matrices
M_NaturalToXYZ = T_xyz*B_natural(:,1:3);
M_XYZToNatural = inv(M_NaturalToXYZ);
M_NaturalToLMS = T_cones*B_natural(:,1:3);
M_LMSToNatural = inv(M_NaturalToLMS);

%% Get background
%
% Start by defining a smooth spectrum with the chromaticity
% we want.  This is not necessarily within the span of the
% subprimaries we have, but serves as a starting point.
targetBgxyY = [targetBgxy ; 1];
targetBgXYZ = xyYToXYZ(targetBgxyY);
desiredBgSpd = B_natural(:,1:3)*M_XYZToNatural*targetBgXYZ;
desiredBgSpd = sum(halfOnSpd)*desiredBgSpd/sum(desiredBgSpd);
desiredBgXYZ = T_xyz*desiredBgSpd;
desiredBgxyY = XYZToxyY(desiredBgXYZ);

%% Search for desired background
%
% Make a guess at reasonable primary values to start at.
startingBgPrimaries = SpdToPrimary(subprimaryCalObjs{1},desiredBgSpd);
startingBgSpd = PrimaryToSpd(subprimaryCalObjs{1},startingBgPrimaries);
startingBgXYZ = T_xyz*startingBgSpd;
startingBgxyY = XYZToxyY(startingBgXYZ);

% You can play with targetLambda a bit by hand to determine tradeoff between
% smoothness and obtaining desired chromaticity.
primaryHeadRoom = 0;
targetLambda = 5;
[bgPrimariesIncr] = ReceptorIsolateSpectral(T_xyz,desiredBgXYZ,subprimaryCalObjs{1}.get('P_device'),startingBgPrimaries,startingBgPrimaries, ...
    primaryHeadRoom,B_natural,projectIndices,targetLambda,subprimaryCalObjs{1}.get('P_ambient'),'EXCITATIONS',true);
bgPrimaries = startingBgPrimaries + bgPrimariesIncr;
bgSettings = PrimaryToSettings(subprimaryCalObjs{1},bgPrimaries);
bgPrimariesQuantized = SettingsToPrimary(subprimaryCalObjs{1},bgSettings);
bgSpd = PrimaryToSpd(subprimaryCalObjs{1},bgPrimaries);
bgLMS = T_cones*bgSpd;
bgXYZ = T_xyz*bgSpd;
bgxyY = XYZToxyY(bgXYZ);
figure; clf; hold on
plot(wls,startingBgSpd,'r','LineWidth',4);
plot(wls,bgSpd,'k','LineWidth',2);
fprintf('Desired  background x,y = %0.3f,%0.3f\n',desiredBgxyY(1),desiredBgxyY(2));
fprintf('Starting background x,y = %0.3f,%0.3f\n',startingBgxyY(1),startingBgxyY(2));
fprintf('Obtained background x,y = %0.3f,%0.3f\n',bgxyY(1),bgxyY(2));
fprintf('Mean value of background primaries: %0.2f\n',mean(bgPrimaries));

end
