% MakeHistogramFlipTime.
%
% This to make a histogram plot using flip time data.
%
% History:
%    08/19/22    smo    - Wrote it.
%    08/26/22    smo    - File name changed.

%% Initialize.
clear; close all;

%% Load the flip time data.
fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/CheckScreenFlipTime';
fileName = '0817_FlipTime_EarlyFrame_0.5.mat';
load(fullfile(fileDir,fileName));

%% Set some parameters.
frequency = 120;
ifi = 1/frequency;

whichRange = 7;

desiredTimeSet = [0.0664 0.4980];
if (whichRange == 6)
    desiredTime = desiredTimeSet(2);
else
    desiredTime = desiredTimeSet(1);
end

%% Plot it.
PLOTREFERENCE = true;

figure; hold on;

% Measured time.
histogram(flipTimeInterval(whichRange,:),'LineWidth',3,'facecolor','b','edgecolor','b');
xlabel('Flip time (sec)', 'fontsize', 15);
ylabel('Counts', 'fontsize', 15);

if (PLOTREFERENCE)
    % Desired time.
    numMaxCount = 120;
    plot([desiredTime desiredTime],[0 numMaxCount],'r-','linewidth',3);
    
    % Desired time +/- one ifi.
    ratioIFI = 1;
    axisIFI = 2.5;
    plot([desiredTime-ratioIFI*ifi desiredTime-ratioIFI*ifi], [0 numMaxCount], 'g--', 'linewidth', 3);
    plot([desiredTime+ratioIFI*ifi desiredTime+ratioIFI*ifi], [0 numMaxCount], 'g--', 'linewidth', 3);
    xlim([desiredTime-axisIFI*ifi desiredTime+axisIFI*ifi]);
    legend('Measured', 'Desired', sprintf('Desired +/- %.1f ifi \n',ratioIFI), 'FontSize', 15);
end
