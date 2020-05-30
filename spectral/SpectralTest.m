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
theIsolatingNaturalSpd = 1.2*sum(halfOnSpd)*theIsolatingNaturalSpd/sum(theIsolatingNaturalSpd);

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
basisType = 'fourier';
switch (basisType)
    case 'cieday'
        load B_cieday
        theNaturalBasis = SplineSpd(S_cieday,B_cieday,S);
    case 'fourier'
        theNaturalBasis = MakeFourierBasis(S,7);
    otherwise
        error('Unknown basis set specified');
end

% Transform
M_NaturalPrimariesToCones = Tcones*theNaturalBasis(:,1:size(Tcones,1));
M_ConesToNaturalPrimaries = inv(M_NaturalPrimariesToCones);

% Get background within basis
theBgNaturalPrimaries = theNaturalBasis\theIsolatingNaturalSpd;
theBgNaturalSpd = theNaturalBasis*theBgNaturalPrimaries;
theBgNaturalLmsTarget = Tcones*theBgNaturalSpd;
figure; clf; hold on
plot(wls,theBgNaturalSpd,'b:','LineWidth',2);

%% Approximation to desired background
IDENTITY = false;
if (IDENTITY)
    B_DevicePrimary = eye(size(cal.computed.pr650M,1));
    theBgDevicePrimariesApprox = theBgNaturalSpd;
    initialDevicePrimariesApprox = theIsolatingNaturalSpd;
else
    B_DevicePrimary = cal.computed.pr650M;
    theBgDevicePrimariesApprox = OLSpdToPrimary(cal,theBgNaturalSpd);
    initialDevicePrimariesApprox = theBgDevicePrimariesApprox;
end
fprintf('Mean value of background primaries: %0.2f\n',mean(theBgDevicePrimariesApprox));

%% Search for desired background
ambientSpd = cal.computed.pr650MeanDark;
primaryHeadRoom = 0;
targetBasis = theNaturalBasis;
targetLambda = 10;
[backgroundDevicePrimariesIncr] = ReceptorIsolateSpectral(Tcones,theBgNaturalLmsTarget',B_DevicePrimary,theBgDevicePrimariesApprox,initialDevicePrimariesApprox, ...
    primaryHeadRoom,targetBasis,projectIndices,targetLambda,ambientSpd,'EXCITATIONS',true);
theBgDevicePrimaries = theBgDevicePrimariesApprox + backgroundDevicePrimariesIncr;
initialDevicePrimaries = theBgDevicePrimaries;
theBgDeviceSpd = OLPrimaryToSpd(cal,theBgDevicePrimaries);
theBgDeviceLms = Tcones*theBgDeviceSpd;
fprintf('Desired/obtained background excitations\n');
for rr = 1:length(theBgNaturalLmsTarget)
    fprintf('\tReceptor %d (desired/obtained): %0.2f, %0.2f\n',rr,theBgNaturalLmsTarget(rr),theBgDeviceLms(rr));
end

% Define desired cone contrast
theLmsContrast = [1 -1 0]'/10;
targetLambda = 10;
[isolatingModulationDevicePrimaries] = ReceptorIsolateSpectral(Tcones,theLmsContrast',B_DevicePrimary,theBgDevicePrimaries,initialDevicePrimaries, ...
    primaryHeadRoom,targetBasis,projectIndices,targetLambda,ambientSpd);
theBgDeviceSpd = B_DevicePrimary*theBgDevicePrimaries;
isolatingDevicePrimariesUpper = isolatingModulationDevicePrimaries + theBgDevicePrimaries;
isolatingDevicePrimariesLower = -isolatingModulationDevicePrimaries + theBgDevicePrimaries;
theIsolatingDeviceSpd = OLPrimaryToSpd(cal,isolatingDevicePrimariesUpper);
theBgDeviceLMS = Tcones*theBgDeviceSpd;
theIsolatingDeviceLMS = Tcones*theIsolatingDeviceSpd;
theIsolatingContrast = ExcitationsToContrast(theIsolatingDeviceLMS,theBgDeviceLMS);
fprintf('Desired/obtained contrasts\n');
for rr = 1:length(theLmsContrast)
    fprintf('\tReceptor %d (desired/obtained): %0.2f, %0.2f\n',rr,theLmsContrast(rr),theIsolatingContrast(rr));
end
fprintf('Min/max primaries upper: %0.4f, %0.4f, lower: %0.4f, %0.4f\n', ...
    min(isolatingDevicePrimariesUpper), max(isolatingDevicePrimariesUpper), ...
    min(isolatingDevicePrimariesLower),  max(isolatingDevicePrimariesLower));

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



