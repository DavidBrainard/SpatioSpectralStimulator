% HistogramFlipTime.
%
% This to make a histogram plot using flip time data.
%
% History:
%    08/19/22    smo    - Wrote it.

%% Initialize.
clear; close all;

%% Load the flip time data.
fileDir = '/Users/seminoh/Aguirre-Brainard Lab Dropbox/Semin Oh/SACC_materials/CheckScreenFlipTime';
fileName = '0817_FlipTime_EarlyFrame_0.5.mat';
load(fullfile(fileDir,fileName));

%% Set some parameters.
frequency = 120;
ifi = 1/frequency;

desiredTimeSet = [0.0664 0.4980];
desiredTime = desiredTimeSet(1);

whichRange = 2;

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
    ratioIFI = 0.5;
    plot([desiredTime-ratioIFI*ifi desiredTime-ratioIFI*ifi], [0 numMaxCount], 'g--', 'linewidth', 3);
    plot([desiredTime+ratioIFI*ifi desiredTime+ratioIFI*ifi], [0 numMaxCount], 'g--', 'linewidth', 3);
    xlim([desiredTime-ratioIFI*ifi desiredTime+ratioIFI*ifi]);
    legend('Measured', 'Desired', sprintf('Desired +/- %.1f ifi \n',ratioIFI), 'FontSize', 15);
end
