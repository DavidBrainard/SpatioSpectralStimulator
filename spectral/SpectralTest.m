% SpectralTest
% 
% Start exploring spectral fits with LEDs
% 
% 4/22/2020  Started on it

%% Clear
clear; close all;

%% Get basis
S = [380 1 401];
wls = SToWls(S);
ledBasis = sssGetSpectralBasis(S);

%% Fit a spectrum with the basis
theTemp = 6500;
theSpd = GenerateBlackBody(theTemp,wls);
weights = lsqnonneg(ledBasis,theSpd);
theFitSpd = ledBasis*weights;
figure; clf; hold on
plot(wls,theSpd,'k','LineWidth',3);
plot(wls,theFitSpd,'r','LineWidth',2);

