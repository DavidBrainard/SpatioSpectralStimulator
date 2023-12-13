%% SquareAndSine
%
% Explicitly show how to decompose a square wave into
% a discrete sum of sinusoids (Fourier series).  Nothing
% deep here, but want to make sure we understand how
% the contrast of fit components comes out when we do
% this.
%
% This link is helpful for the formula.
%  https://mathworld.wolfram.com/FourierSeriesSquareWave.html

%% Initialize
clear; close all;

%% Define square wave parameters
nSamples = 1000;
maxSpace = 100;
x = linspace(0,maxSpace,nSamples);
sf = 4;
period = maxSpace/sf;
contrast = 0.9;
mean = 3;

%% Build and plot square wave
sw = square(2*pi*x/period);
sw = sw*contrast*mean + mean;
figure; clf; hold on;
plot(x,sw,'k','LineWidth',4);
actualContrast = (max(sw)-min(sw))/(max(sw)+min(sw));
fprintf('Desired contrast = %0.4f, actual = %0.4f\n',contrast,actualContrast);

%% Camera MTF simulation
% 
% If you set this to true, then there is a rolloff of the
% sinusoidal component contrasts with spatial frequency.
% In this case, the (max-min)/(max+min) estimate of contrast
% from the nominal square wave will increase relative to the sinusoidal
% contrast of the fundamental, which itself doesn't change.
% This illustrates the effect we were talking about today.
% You see it by comparing the last number printet out at the
% end with the numbers in the line before.
%
% Note that you need nSines to be large for this to work, otherwise
% the square wave approx rings too much for a simple extraction
% of contrast to work well.
SIMMTF = false;
simMTFExp = 0.2;

%% Build up fourier series
nSines = 1000;
swApprox = mean;
sineComponentIndex = 0;
for ss = 1:2:2*(nSines-1)+1
    sineComponentIndex = sineComponentIndex + 1;

    if (SIMMTF)
        MTFFACTOR(sineComponentIndex) = (1/ss)^simMTFExp;
    else
        MTFFACTOR(sineComponentIndex) = 1;
    end
     
    sineComponent{sineComponentIndex} = contrast*mean*(4/pi)*(1/ss)*sin(ss*2*pi*x/period);
    swApprox = swApprox + MTFFACTOR(sineComponentIndex)*sineComponent{sineComponentIndex};

    sineComponentPredContrast(sineComponentIndex) = contrast*(4/pi)*(1/ss);
    sineComponentActualContrast(sineComponentIndex) = (max(sineComponent{sineComponentIndex})-min(sineComponent{sineComponentIndex}))/(max(sineComponent{sineComponentIndex})+min(sineComponent{sineComponentIndex})+2*mean);
end
plot(x,swApprox,'r','LineWidth',2);
plot(x,sineComponent{1}+mean,'g');
xlabel('Position');
ylabel('Amplitude');
title(sprintf('%d c/image, contrast %0.2f',sf,contrast));

%% Plot predicted versus "measured" component contrasts
figure; clf; hold on;
plot(sineComponentPredContrast,sineComponentActualContrast,'ro','MarkerFaceColor','r','MarkerSize',8);
plot([0 1.5],[0 1.5],'k:')
axis('square');
xlabel('Predicted Component Contrast Around Mean');
ylabel('Actual Component Contrast Around Mean');

%% Contrast of fundamental sine 
fprintf('Contrast of fundamental sinusoid is %0.4f. 4/pi*contrast is %0.4f\n',sineComponentActualContrast(1),sineComponentPredContrast(1));

%% Relation between measured square wave contrast and fundamental sine contrast
swApproxContrast = (max(swApprox)-min(swApprox))/(max(swApprox)+min(swApprox));
fprintf('Square wave approx contrast %0.4f, times 4/pi %0.4f\n',swApproxContrast,4/pi*swApproxContrast)