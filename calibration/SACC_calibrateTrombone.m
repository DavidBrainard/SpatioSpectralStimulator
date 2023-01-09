% SACC_calibrateTrombone
%
% This contains the contrast measurement data over Tromobone position.
% Based on the data, you can search an optimal Trombone position for any
% given lens diopters.

% History:
%    03/25/22  dhb, smo     - Measured the data and plotted the graph.
%    07/19/22  dhb, smo     - Added data using different lense diopters.
%    07/29/22  smo, rsh     - Added data using different lenese within the
%                             range of -6 to +2.
%    08/01/22  smo          - Added to choose multiple data at once. Also,
%                             now it calculates an optimal Trombone
%                             position for any given lens diopters using
%                             interpolation.
%    08/02/22  smo          - Changed the fitting method from interpolation
%                             to linear. And file name changed.
%    10/20/22  smo          - New calibration data added for the lens
%                             diopters of +2 to +6.

%% Initialize.
clear; close all;

%% Raw measurement data.
%
% Trombone position was recorded in mm unit. For the data (whichData =
% 0-10), the position of Trombone is the number you can get from its
% software. However, that number is changing everytime rebooting the
% computer which does not represent the physical position. So, the data
% (whichData = 11-18) recorded the number on the measuring tape attached to
% the side of Trombone, which is an actual position of the Trombone.
%
% Contrast was measured using the real-time camera image.
%
% Note that we are planning to calibrate for the lens diopeters from -6.0
% to +2.0 for SACC project.
%
% Set the data range of interest.
whichDataRange = [11:22];
VERBOSE = false;

% Make a loop here to collect the multiple lens data.
for ii = 1:length(whichDataRange)
    whichData = whichDataRange(ii);
    
    % All raw data stored here.
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
            
            % From here, data was collected on the day of 07/29/22.
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
            %         positionTrombone = [130.5 129.5 127.5 126 125 122.5 121];
            %         contrast         = [0.909 0.927 0.938 0.939 0.938 0.924 0.883];
            %
            % Data for -4.0 has been corrected as follows. The first
            % measurement somehow gave the lower contrast than expected.
            positionTrombone = [130.5 129.5 128 127 126 125 124 122.5 121.5 120.5];
            contrast         = [0.901 0.923 0.933 0.946 0.955 0.944 0.943 0.933 0.92 0.894];
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
        case 17
            lensDiopters = +1.0;
            positionTrombone = [165 164 162.5 161.5 160.5 160 159 157.5 156 155];
            contrast         = [0.901 0.942 0.95 0.959 0.959 0.959 0.958 0.95 0.941 0.872];
        case 18
            lensDiopters = +2.0;
            positionTrombone = [170.5 169.5 168.5 168 166 165 164 162 161];
            contrast         = [0.908 0.938 0.953 0.953 0.961 0.961 0.953 0.937 0.896];
            
            % This is the newly obtained data to validate the measurement
            % (as of 10/20/22). We repeated the meausurement of lens
            % diopters of +2 to see if we can get the same results as the
            % last time and it did!
        case 23
            lensDiopters = +2.0;
            positionTrombone = [170.5 169.5 168.5 168 166 165 164 162 161];
            contrast         = [0.908 0.935 0.943 0.942 0.956 0.956 0.949 0.922 0.912];
        case 19
            lensDiopters = +3.0;
            positionTrombone = [176 175 174 173 172 171 170 169 168 167 166 165];
            contrast         = [0.897 0.938 0.937 0.945 0.950 0.950 0.952 0.957 0.944 0.935 0.900 0.841];
        case 20
            lensDiopters = +4.0;
            positionTrombone = [182 181 180 179 178 177 176 175 174 173 172 171];
            contrast         = [0.881 0.914 0.936 0.943 0.951 0.951 0.958 0.952 0.944 0.936 0.914 0.879];
        case 21
            lensDiopters = +5.0;
            positionTrombone = [186 185 184 183 182 181 180 179 178 177 176 175];
            contrast         = [0.894 0.923 0.942 0.949 0.949 0.953 0.95 0.942 0.936 0.934 0.905 0.831];
        case 22
            lensDiopters = +6.0;
            positionTrombone = [191 190 189 188 187 186 185 184 183 182 181 180];
            contrast         = [0.889 0.910 0.931 0.939 0.947 0.954 0.945 0.944 0.935 0.921 0.914 0.861];
        otherwise
    end
    
    % There are some lenses having multiple highest points. Here we choose a
    % smaller value in position (nm) which should be fine for now, but we
    % may elaborate it later on.
    lensDioptersSet(ii) = lensDiopters;
    contrastMax(ii) = max(contrast);
    positionContrastMax(ii) = min(positionTrombone(find(contrast(:) == contrastMax(ii))));
    
    %% Plot it.
    if (VERBOSE)
        figure; clf; hold on;
        plot(positionTrombone, contrast, 'ko--', ...
            'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', zeros(1,3), 'MarkerSize', 8);
        plot(positionContrastMax(ii), contrastMax(ii), 'o', ...
            'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', zeros(1,3), 'MarkerSize', 8);
        xlabel('Position Trombone (mm)','FontSize',15);
        ylabel('Contrast','FontSize',15);
        xticks([round(min(positionTrombone)) : 1 : round(max(positionTrombone))]);
        ylim([0.855 0.965]);
        legend('All data','Max Contrast','location','southeast','FontSize',15);
        title(append('Lens diopters: ',num2str(lensDiopters),' / Max contrast = ',num2str(contrastMax(ii))),'FontSize',15);
    end
    
    clear contrast; clear positionTrombone;
end

%% Search an optimal Trombone position for a random lens diopters.
%
% Type your lens diopter of interest.
PROMPT = true;

if (PROMPT)
    while 1
    promptMessage = '    Type subject lens diopter of dominant eye: ';
    subjectLensDiopters = input(promptMessage);
    
    if ~isempty(subjectLensDiopters)
        break
    end
    
    disp('Lens diopter should be within from -6 to +2.');
    end
else
    % Enter manually.
    subjectLensDiopters = -5.25;
end

% Set the target (initial) position.
targetLensDiopters = -(subjectLensDiopters + 0.5);

% Here make a linear fitting for the data.
p = polyfit(lensDioptersSet, positionContrastMax, 1);
positionFitted = polyval(p ,lensDioptersSet);

% Search an optimal Trombone position here from the above fitting results.
positionEstimated = polyval(p, targetLensDiopters);

% Print out the results.
fprintf('\t * Subject lens diopter is (%.2f) \n', subjectLensDiopters);
fprintf('\t * Initial Trombone position at target lens diopter (%.2f) is %.0f (mm) \n', ...
    targetLensDiopters, positionEstimated);

% Plot it.
PLOTCALIBRATIONRESULTS = true;

if (PLOTCALIBRATIONRESULTS)
    figure; clf; hold on;
    plot(lensDioptersSet, positionContrastMax, 'k.', 'markersize', 20);
    plot(lensDioptersSet, positionFitted, 'r--');
    plot(targetLensDiopters, positionEstimated, ...
        'ro', 'markersize', 13, 'markerfacecolor', 'r', 'markeredgecolor', 'k');
    xlabel('Testing lens diopters','FontSize',15);
    ylabel('Trombone position (mm)','FontSize',15);
    title('Initial Trombone position over lense diopers','FontSize',15);
    legend('Measure','Fit', 'Estimation', 'location', 'SouthEast', 'FontSize', 12);
end
