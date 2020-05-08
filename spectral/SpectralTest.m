% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Get basis
% S = [380 1 401];
% wls = SToWls(S);
% ledBasis = sssGetSpectralBasis(S);

%% Generate a OneLight cal file
S = [400 2 151];
cal = sssSpoofOneLightCal('S',S, ....
    'plotBasis',true,...
    'gaussianPeakWls',[437 485 540 562 585 618 652], ...
    'gaussianFWHM',25);
S = cal.describe.S;
wls = SToWls(S);

%% Cone fundamentals
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
Tcones = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);

%% Get half on spectrum
halfOnPrimaries = 0.5*ones(cal.describe.numWavelengthBands,1);
halfOnSpd = OLPrimaryToSpd(cal,halfOnPrimaries);
figure; hold on;
plot(wls,halfOnSpd,'r','LineWidth',3);

%% Make sure gamma correction behaves well
halfOnSettings = OLPrimaryToSettings(cal,halfOnPrimaries);
if (max(abs(halfOnPrimaries-halfOnSettings)) > 1e-8)
    error('Spoofed cal file does not have linear gamma');
end
halfOnPrimariesChk = OLSettingsToPrimary(cal,halfOnSettings);
if (max(abs(halfOnPrimaries-halfOnPrimariesChk)) > 1e-8)
    error('Spoofed cal file does not have linear inverse gamma');
end

%% Use OL machinery to get primaries from spectrum
halfOnPrimariesChk = OLSpdToPrimary(cal,halfOnSpd);
halfOnSpdChk = OLPrimaryToSpd(cal,halfOnPrimariesChk);
plot(wls,halfOnSpdChk,'k','LineWidth',1);

%% Fit a spectrum with the basis
%
% Generate a black body radiator and put it into right general scale
theTemp = 6500;
theSpd = GenerateBlackBody(theTemp,wls);
theSpd = 0.6*sum(halfOnSpd)*theSpd/sum(theSpd);

% Use OL machinery to get primaries and then the spd
thePrimaries = OLSpdToPrimary(cal,theSpd,'lambda',0.01);
theFitSpd = OLPrimaryToSpd(cal,thePrimaries);

% Use non-neg least squares for comparison
weights = lsqnonneg(cal.computed.pr650M,theSpd);
theFitSpdRegress = cal.computed.pr650M*weights;

% Plot
figure; clf; hold on
plot(wls,theSpd,'k','LineWidth',3);
plot(wls,theFitSpd,'r','LineWidth',2);
plot(wls,theFitSpdRegress,'g','LineWidth',1);

%% Generate cone isolating spectral modulations with a simple basis
load B_cieday
theBasis = SplineSpd(S_cieday,B_cieday,S);
M_BasisToCones = Tcones*theBasis;
M_ConesToBasis = inv(M_BasisToCones);

% Get background within basis
theBgPrimaries = lsqnonneg(theBasis,theSpd);
theBgSpd = theBasis*theBgPrimaries;
theBgLms = Tcones*theBgSpd;
plot(wls,theBgSpd,'b:','LineWidth',2);

% Define desired cone contrast
theLmsContrast = [1 -1 0]';
theLmsInc = theLmsContrast.*theBgLms;
thePrimariesDir = M_ConesToBasis*theLmsInc;
thePrimariesInc = 
theLms = theLmsInc+theBgLms;
thePrimaries = M_ConesToBasis*theLms;
theSpd = theBasis*thePrimaries;
plot(wls,theSpd,'b','LineWidth',2);

% modulationBasis = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
%     theBasis, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
%     primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

