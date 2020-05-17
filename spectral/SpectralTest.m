% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Generate a OneLight cal file
S = [380 2 201];
cal = sssSpoofOneLightCal('S',S, ....
    'plotBasis',true,...
    'gaussianPeakWls',[437 485 540 562 585 618 652], ...
    'gaussianFWHM',25);
S = cal.describe.S;
wls = SToWls(S);
lowProjectWl = 420;
highProjectWl = 680;
projectIndices = find(wls > lowProjectWl & wls < highProjectWl);

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
theIsolatingNaturalSpd = GenerateBlackBody(theTemp,wls);
theIsolatingNaturalSpd = 0.6*sum(halfOnSpd)*theIsolatingNaturalSpd/sum(theIsolatingNaturalSpd);

% Use OL machinery to get primaries and then the spd
theIsolatingNaturalPrimaries = OLSpdToPrimary(cal,theIsolatingNaturalSpd,'lambda',0.01);
theFitSpd = OLPrimaryToSpd(cal,theIsolatingNaturalPrimaries);

% Use non-neg least squares for comparison
weights = lsqnonneg(cal.computed.pr650M,theIsolatingNaturalSpd);
theFitSpdRegress = cal.computed.pr650M*weights;

% Plot
figure; clf; hold on
plot(wls,theIsolatingNaturalSpd,'k','LineWidth',3);
plot(wls,theFitSpd,'r','LineWidth',2);
plot(wls,theFitSpdRegress,'g','LineWidth',1);

%% Generate cone isolating spectral modulations with a simple basis
load B_cieday
theNaturalBasis = SplineSpd(S_cieday,B_cieday,S);
M_NaturalPrimariesToCones = Tcones*theNaturalBasis;
M_ConesToNaturalPrimaries = inv(M_NaturalPrimariesToCones);

% Get background within basis
theBgNaturalPrimaries = theNaturalBasis\theIsolatingNaturalSpd;
theBgNaturalSpd = theNaturalBasis*theBgNaturalPrimaries;
theBgNaturalLms = Tcones*theBgNaturalSpd;
figure; clf; hold on
plot(wls,theBgNaturalSpd,'b:','LineWidth',2);

% Define desired cone contrast
theLmsContrast = [1 -1 0]'/10;
theLmsInc = theLmsContrast.*theBgNaturalLms;
theLms = theLmsInc+theBgNaturalLms;
theIsolatingNaturalPrimaries = M_ConesToNaturalPrimaries*theLms;
theIsolatingNaturalSpd = theNaturalBasis*theIsolatingNaturalPrimaries;
plot(wls,theIsolatingNaturalSpd,'b','LineWidth',2);

%% Now we want to find an in-gamut modulation 
IDENTITY = false;
if (IDENTITY)
    B_DevicePrimary = eye(size(cal.computed.pr650M,1));
    theBgDevicePrimaries = theBgNaturalSpd;
    initialDevicePrimaries = theIsolatingNaturalSpd;
else
    B_DevicePrimary = cal.computed.pr650M;
    theBgDevicePrimaries = OLSpdToPrimary(cal,theBgNaturalSpd);
    initialDevicePrimaries = theBgDevicePrimaries;
end

whichReceptorsToIsolate = [1 2];
whichReceptorsToIgnore = [];
whichReceptorsToZero = setdiff(1:size(Tcones,1),[whichReceptorsToIsolate whichReceptorsToIgnore]);
desiredContrast = [1 -1]/10;

ambientSpd = cal.computed.pr650MeanDark;
maxPowerDiff = Inf;
targetBasis = theNaturalBasis;
targetLambda = 2;
primaryHeadRoom = 0;

[isolatingModulationDevicePrimaries] = ReceptorIsolateSpectral(Tcones,whichReceptorsToIsolate, ...
    whichReceptorsToIgnore,B_DevicePrimary,theBgDevicePrimaries,initialDevicePrimaries, ...
    primaryHeadRoom,maxPowerDiff,targetBasis,projectIndices,targetLambda,desiredContrast,ambientSpd);
theBgDeviceSpd = B_DevicePrimary*theBgDevicePrimaries;
isolatingDeviceUpperPrimaries = isolatingModulationDevicePrimaries + theBgDevicePrimaries;
isolatingDevicePrimariesLower = -isolatingModulationDevicePrimaries + theBgDevicePrimaries;
theIsolatingDeviceSpd = B_DevicePrimary*isolatingDeviceUpperPrimaries;
theBgDeviceLMS = Tcones*theBgDeviceSpd;
theIsolatingDeviceLMS = Tcones*theIsolatingDeviceSpd;
theIsolatingContrast = ExcitationsToContrast(theIsolatingDeviceLMS,theBgDeviceLMS);
fprintf('Desired/obtained contrasts\n');
for ii = 1:length(whichReceptorsToIsolate)
    rr = whichReceptorsToIsolate(ii);
    fprintf('\tReceptor %d (desired/obtained): %0.2f, %0.2f\n',rr,desiredContrast(rr),theIsolatingContrast(rr));
end
for ii = 1:length(whichReceptorsToZero)
    rr = whichReceptorsToZero(ii);
    fprintf('\tReceptor %d (desired/obtained): %0.2f, %0.2f\n',rr,0,theIsolatingContrast(rr));
end
for ii = 1:length(whichReceptorsToIgnore)
    rr = whichReceptorsToIgnore(ii);
    fprintf('\tReceptor %d,ignored (obtained): %0.2f\n',rr,theIsolatingContrast(rr));
end

% How close are spectra to subspace defined by basis?
theBgNaturalApproxSpd = targetBasis*(targetBasis(projectIndices,:)\theBgDeviceSpd(projectIndices));
theIsolatingNaturalApproxSpd = targetBasis*(targetBasis(projectIndices,:)\theIsolatingDeviceSpd(projectIndices));

figure; clf; 
subplot(1,2,1); hold on
plot(wls,theBgDeviceSpd,'b','LineWidth',2);
plot(wls,theBgNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),theBgDeviceSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theBgNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
ylim([0 2]);
subplot(1,2,2); hold on
plot(wls,theIsolatingDeviceSpd,'b','LineWidth',2);
plot(wls,theIsolatingNaturalApproxSpd,'r:','LineWidth',1);
plot(wls(projectIndices),theIsolatingDeviceSpd(projectIndices),'b','LineWidth',4);
plot(wls(projectIndices),theIsolatingNaturalApproxSpd(projectIndices),'r:','LineWidth',3);
ylim([0 2]);



