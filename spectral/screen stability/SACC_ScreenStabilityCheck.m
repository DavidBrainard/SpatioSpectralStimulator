% SACC_ScreenStabilityCheck
%
% This tests screen stability for SACC project. It basically measures
% the plain screen repeatedely and 
% 
% History:
%    11/24/2021 smo   Pulled out this part from the old code. It has been
%                     cleaned up using our SACC measurement functions.

%% Set parameters here.
% 
% This code basically measures the spectrum of the plain screen
% automatically in a specific time interval.

% Set measurement time length and interval. The time is set in minute unit,
% which will be converted into second unit later in this code.
totalMeasurementTime_min = 0.2;
timeDelayBeforeMeasurement_min = 0.1;
timeDelayBeforeMeasurement_sec = timeDelayBeforeMeasurement_min * 60;

% Measurement range.
S = [380 2 201];

% Verbose.
verbose = true;

%% Make a bit of time delay before the measurement starts. 
% You can go out before it's done.
timeDelayGoOut_sec = 3;
for tt = 1:timeDelayGoOut_sec
     pause(1) % Pause for 1 seconds (unit)
end

%% Make screen and spectroradiometer ready.
%
% Open the plain screen. Simply set it as white here. It won't change
% during the whole measurements.
screenSettings = [1 1 1];
OpenPlainScreen(screenSettings,'projectorMode',true,'verbose',verbose);

% Connect to the spectroradiometer. We will use PR670 here.
OpenSpectroradiometer;

%% Measurements.
% 
% The number of measurements are calculated based on the total measurement
% time and time delay interval which were set from the above.
%
% +1 is the measurement at the cold state which is right after the screen turned on.
nMeasurments = (totalMeasurementTime_min / timeDelayBeforeMeasurement_min) + 1; 
allSpdMeasured = zeros(S(3),nMeasurments);

% Measure it.
for i=1:nMeasurments
    
    % Measurement happens here.
    allSpdMeasured(:,i) = MeasureSpectroradiometer;
    
    % Make a delay before next measurement.
    for tt = 1:timeDelayBeforeMeasurement_sec
        % Pause for 1 second in unit.
        pause(1);
    end
end

%% XYZ calculations
%
% Match the wavelength range.
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos, 683*T_xyzJuddVos, S);

% Calculate XYZ values.
XYZ = 683 * T_xyz * allSpdMeasured;
xyY = XYZToxyY(XYZ);

% Calculate color gamut too.
colorGamut = XYZToxyY(T_xyz);
colorGamut(:,end+1) = colorGamut(:,1);

%% Plot the data.
if (verbose)
   % Spds.
   figure; clf;
   plot(allSpdMeasured);
   xlabel('Wavelenth (nm)');
   ylabel('Spectral power distribution');

   % Luminance.
   figure; clf;
   subplot(1,2,1);
   measurementTime = linspace(0, totalMeasurementTime_min, nMeasurments);
   plot(measurementTime, XYZ(3,:),'r*--');
   xlabel('Measurement time (min)');
   ylabel('Luminance (cd/m2)');
   ylim([max(XYZ(3,:))*0.8 max(XYZ(3,:))*1.2]);
   legend('Measurements');
   
   % xy coordiantes.
   subplot(1,2,2); hold on;
   plot(xyY(1,:), xyY(2,:), 'r*');
   plot(colorGamut(1,:),colorGamut(2,:),'k-');
   xlabel('CIE x');
   ylabel('CIE y');
   legend('Measurements','Color Gamut');
end

%% Save the data.
%
% Close the screen and projector.
CloseScreen;
CloseSpectroradiometer;

% Save data with the name containing dayTimestr.
if (ispref('SpatioSpectralStimulator','CheckDataFolder'))
    testFiledir = getpref('SpatioSpectralStimulator','CheckDataFolder');
    dayTimestr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
    testFilename = fullfile(testFiledir,sprintf('stabilityCheck_%s',dayTimestr));
    save(testFilename,'allSpdMeasured','XYZ','xyY','colorGamut', ...
                      'measurementTime','S');
end
 