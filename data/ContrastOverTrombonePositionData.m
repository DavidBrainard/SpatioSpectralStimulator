% ContrastOverTrombonePositionData
%
% This is the measurement of contrast over the Tromobone position by
% controlling it by its own software. 

% History:
%    03/25/22  dhb, smo     - Measured the data and plotted the graph.
%    07/19/22  dhb, smo     - Added data using different lenses.

%% Measurement data.
%
% Position is in mm unit in the Trombone controlling software, and the
% contrast was measured based on the real-time camera image.
%
% Note that we are planning to calibrate for the lens diopeters from -6.0
% to +2.0 for SACC project.
whichData = 1;

switch whichData
    case 0
        lensDiopeters = 0;
        positionTrombone = [2.799 2.299 1.772 1.522 1.272 1.022 0.772 0.522 0.272 0.022 -0.530 -1.449 -2.193 -2.665 ];
        contrast         = [0.931 0.946 0.954 0.954 0.954 0.954 0.955 0.954 0.954 0.954  0.954  0.954  0.947  0.939 ];
    % Following data no. 1 to 6 were collected using different lens
    % diopeters on the date of 07/19/22.
    case 1
        lensDiopeters = -3;
        positionTrombone = [-15.17 -14.64 -14.06 -13.54 -12.99 -12.52 -12.05 -11.52 -11.01 -10.54 -9.99 -9.48];
        contrast         = [ 0.888  0.871  0.901  0.929  0.938  0.947  0.948  0.948  0.948  0.949 0.931 0.923];
    case 2
        lensDiopeters = -4.5;
        positionTrombone = [-19.66 -19.54 -19.04 -18.04 -17.50];
        contrast         = [ 0.946  0.937  0.947  0.935  0.923];
    case 3
        lensDiopeters = +2.0;
        positionTrombone = [-0.47  0.54  1.00  1.50  2.06  2.50  3.04];
        contrast         = [0.920 0.947 0.947 0.963 0.946 0.947 0.939];
    case 4
        lensDiopeters = -6.0;
        positionTrombone = [-17.89 -16.80];
        contrast         = [ 0.815  0.718];
    case 5
        lensDiopeters = -5.5;
        positionTrombone = [-19.79 -18.82];
        contrast         = [ 0.917  0.824];
    case 6
        lensDiopeters = -5.0;
        positionTrombone = [-19.53 -18.55];
        contrast         = [ 0.930  0.911];
    otherwise
end
        
contrastMax = max(contrast);
positionContrastMax = positionTrombone(find(contrast(:) == contrastMax));

%% Plot it.
figure; clf; hold on;
plot(positionTrombone, contrast, 'ko--', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', zeros(1,3), 'MarkerSize', 8);
plot(positionContrastMax, contrastMax, 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', zeros(1,3), 'MarkerSize', 8);
xlabel('Position Trombone (mm)','FontSize',13);
ylabel('Contrast','FontSize',13);
ylim([0.871 0.965]);
legend('All data','Max Contrast','location','southeast');
title(append('Lens diopters: ',num2str(lensDiopeters),' / Max contrast = ',num2str(contrastMax)));
