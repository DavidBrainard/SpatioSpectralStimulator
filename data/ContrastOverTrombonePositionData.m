% ContrastOverTrombonePositionData
%
% This is the measurement of contrast over the Tromobone position by
% controlling it by its own software. 

% History:
%    03/25/22  dhb, smo     - Measured the data and plotted the graph.
%    07/19/22  dhb, smo     - Added data using different lense diopters.
%    07/29/22  smo, rsh     - Added data using different lenese within the
%                             range of -6 to +2.

%% Measurement data.
%
% Position is in mm unit in the Trombone controlling software, and the
% contrast was measured based on the real-time camera image.
%
% Note that we are planning to calibrate for the lens diopeters from -6.0
% to +2.0 for SACC project.
whichData = 11;

switch whichData
    case 0
        lensDiopters = 0;
        positionTrombone = [2.799 2.299 1.772 1.522 1.272 1.022 0.772 0.522 0.272 0.022 -0.530 -1.449 -2.193 -2.665 ];
        contrast         = [0.931 0.946 0.954 0.954 0.954 0.954 0.955 0.954 0.954 0.954  0.954  0.954  0.947  0.939 ];
    % Following data no. 1 to 6 were collected using different lens
    % diopeters on the date of 07/19/22.
    case 1
        lensDiopters = -3;
        positionTrombone = [-15.17 -14.64 -14.06 -13.54 -12.99 -12.52 -12.05 -11.52 -11.01 -10.54 -9.99 -9.48];
        contrast         = [ 0.888  0.871  0.901  0.929  0.938  0.947  0.948  0.948  0.948  0.949 0.931 0.923];
    case 2
        lensDiopters = -4.5;
        positionTrombone = [-19.66 -19.54 -19.04 -18.04 -17.50];
        contrast         = [ 0.946  0.937  0.947  0.935  0.923];
    case 3
        lensDiopters = +2.0;
        positionTrombone = [-0.47  0.54  1.00  1.50  2.06  2.50  3.04];
        contrast         = [0.920 0.947 0.947 0.963 0.946 0.947 0.939];
    case 4
        lensDiopters = -6.0;
        positionTrombone = [-17.89 -16.80];
        contrast         = [ 0.815  0.718];
    case 5
        lensDiopters = -5.5;
        positionTrombone = [-19.79 -18.82];
        contrast         = [ 0.917  0.824];
    case 6
        lensDiopters = -5.0;
        positionTrombone = [-19.53 -18.55];
        contrast         = [ 0.930  0.911];
        
        %% From here, data was collected on the day of 07/29/22.
        %
        % Note that the position trombone here is the number on the
        % measuring tape that is attached to the Trombone. Therefore, this
        % is different range of unit from the measurements conducted
        % before. The unit for the position is still mm though.
        
    case 11
        lensDiopters = -6.0;
        positionTrombone = [115 114 112.5 111 110 109 107.5 105];
        contrast         = [0.92 0.947 0.947 0.959 0.96 0.946 0.931 0.914];
    case 12
        lensDiopters = -5.0;
        positionTrombone = [124 121.5 120.5 120 119 117.5 115.5 114.5 112.5];
        contrast         = [0.92 0.96 0.96 0.959 0.959 0.958 0.958 0.928 0.898];
    case 13
        lensDiopters = -4.0;
        positionTrombone = [130.5 129.5 127.5 126 125 122.5 121];
        contrast         = [0.909 0.927 0.938 0.939 0.938 0.924 0.883];
    case 14
        lensDiopters = -3.0;
        positionTrombone = [138.5 137 136 135 134 132.5 131 129 128];
        contrast         = [0.91 0.94 0.941 0.949 0.95 0.95 0.9 0.901 0.856];
    case 15
        lensDiopters = -2.0;
        positionTrombone = [146 145 144 143 141 140 139 138 135.5];
        contrast         = [0.88 0.931 0.94 0.95 0.951 0.959 0.949 0.939 0.855];
    case 16
        lensDiopters = -1.0;
        positionTrombone = [152.5 151.5 150.5 149.5 148.5 147 145.5 144 143];
        contrast         = [0.893 0.929 0.946 0.946 0.947 0.955 0.955 0.946 0.927];
    case 1718
        lensDiopters = +1.0;
        positionTrombone = [165 164 162.5 161.5 160.5 160 159 157.5 156 155];
        contrast         = [0.901 0.942 0.95 0.959 0.959 0.959 0.958 0.95 0.941 0.872];
    case 18
        lensDiopters = +2.0;
        positionTrombone = [170.5 169.5 168.5 168 166 165 164 162 161 159.5];
        contrast         = [0.908 0.938 0.953 0.953 0.961 0.961 0.953 0.937 0.896 0.8];
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
ylim([0.86 0.965]);
legend('All data','Max Contrast','location','southeast');
title(append('Lens diopters: ',num2str(lensDiopters),' / Max contrast = ',num2str(contrastMax)));
