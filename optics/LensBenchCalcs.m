%% LensBenchCalcs
% A simple optical bench explorer
%
% Description:
%   The origin of the optical system [0,0,0] corresponds to the apex of the
%   cornea of the eye. The optical axis of the eye runs along the x
%   dimension. The eye is positioned at negative x values. This routine
%   allows one to specify the position of bi-convex lenses centered along
%   the optical axis. You must add the lenses in the left to right
%   direction. That is, each subsequent lens must have a larger (more
%   positive) position.
%
%   This version developed from DEMO_lensBench in gkaModelEye.

%% Housekeeping
clear; close all

%% Plot options
plotOutputRays = false;

%% Eye params
sphericalAmetropia = 0;
mmPerDeg = 0.3;

%% Set the position of DLP chip
%
% DLP distance from the cornea in mm
totalPixels = 1920;

%% Setup lens params
lensDiameterMm = 50;                  % Lens diameters in mm
lensFocalLengthMm = 200;              % Lens focal lengths in mm
specifiedLensFocalLengthMm = 1.01*lensFocalLengthMm;
nLenses = 3;

% Compute lens center and DLP position from focal lengths, and
% fill in arrays
for ll = 1:nLenses
    lensDiametersMm(ll) = lensDiameterMm;
    lensFocalLengthsMm(ll) = specifiedLensFocalLengthMm;
    lensCenters(ll) = lensFocalLengthMm + (ll-1)*2*lensFocalLengthMm;
end
DLPposition = lensFocalLengthMm*nLenses*2;

%% Stops in system
irisStopRadiusMm = 2;
artificialPupilRadiusMm = 1;
artificalPupilPositionMm = 4*lensFocalLengthMm;

% Convert to useful format
for ll = 1:length(lensCenters)
    lensRadiiMm(ll) = lensDiametersMm(ll)/2;
    lensPowersDiopters(ll) = 1000./lensFocalLengthsMm(ll); 
end

%% Initialize the optical system with an eye
% Create an initial optical system in the eyeToCamera (left to right)
% direction with an emmetropic right eye focused at 1.5 meters, with the
% refractive indices for the visible spectrum.
fprintf('Setup base optical system\n');
sceneGeometry = createSceneGeometry('spectralDomain','vis','sphericalAmetropia',sphericalAmetropia);

%% Add an iris stop
opticalSystemStruct = sceneGeometry.refraction.retinaToCamera;
opticalSystemStruct = addStopAfter(opticalSystemStruct, irisStopRadiusMm);
   
%% Add the lenses to the optical system at desired locations
fprintf('Add lenses\n');
for ll = 1:length(lensCenters)
    opticalSystemStruct = addBiconvexLens( opticalSystemStruct, lensCenters(ll), lensPowersDiopters(ll), lensRadiiMm(ll));
end

%% Add artificial pupil before the last lens
targetRow = length(opticalSystemStruct.surfaceLabels)-2;
opticalSystemStruct = addStopAfter(opticalSystemStruct, artificialPupilRadiusMm, artificalPupilPositionMm, targetRow);

%% Reverse the optical system direction
% The optical system is assembled for ray tracing from left-to-right (away
% from the eye), but we want to do ray tracing from right-to-left (towards
% the eye)
opticalSystemStruct = reverseSystemDirection(opticalSystemStruct);

%% Display the optical system
fprintf('Setup plot\n');
figure;
set(gcf,'Position',[100 725 1290 620]);
%plotOpticalSystem('newFigure',true,'surfaceSet',opticalSystemStruct,'addLighting',true,'viewAngle',[0 90]);

%% Extract the opticalSystem matrix
opticalSystem = opticalSystemStruct.opticalSystem;

%% Add some rays
% Send rays from the horizontal bounds of the DLP chip, which will be 11.5
% mm on either side of the optical axis, and have an angle w.r.t. the
% optical axis of +- 12 degrees
numberRaysPerBundle = 400;
rayAnglesDeg = linspace(-12,12,numberRaysPerBundle);
rayHorizStartPosMm = [-10.5,0,10.5];
colors = {'red','green','blue'};

% Loop over the desired rays and display
retinalHorizPositions = zeros(length(rayHorizStartPosMm), length(rayAnglesDeg));
for hh = 1:length(rayHorizStartPosMm)
    fprintf('Trace ray bundle %d of %d\n',hh,length(rayHorizStartPosMm));
    color = colors{hh};
    for aa = 1:length(rayAnglesDeg)
        % Create the ray
        inputRay = quadric.normalizeRay(quadric.anglesToRay([DLPposition;rayHorizStartPosMm(hh);0],-180+rayAnglesDeg(aa),0));
        
        % Trace it
        %
        % For the location in the retina: The outputRay of the ray trace is
        % of the form [p; d], where p is the point of intersection on the
        % retina in [x,y,z] coords [axial,horizontal,vertical]
        [outputRay, rayPath] = rayTraceQuadrics(inputRay, opticalSystem);
        retinalHorizPositions(hh,aa) = outputRay(2,1);
        
        % Add it to the plot
        plotOpticalSystem('newFigure',false,'rayPath',rayPath,'rayColor',color,'viewAngle',[0 90]);
        
        % If the outputRay is nan (that is, the ray missed the eye), then
        % extend the ray as it left the last lens surface
        if any(any(isnan(outputRay))) && plotOutputRays
            surfaces = find(~any(isnan(rayPath)));            
            outputRay = rayTraceQuadrics(inputRay, opticalSystem(surfaces,:));
            outputRay(:,2) = outputRay(:,2) * 5; % Pump up the volume
            plotOpticalSystem('newFigure',false,'outputRay',outputRay,'rayColor',color,'viewAngle',[0 90]);  
        end

    end
end

%% Find retinal extent
retinalLocationsMm = nanmean(retinalHorizPositions,2);
retinalExtentMm = max(retinalLocationsMm)-min(retinalLocationsMm);
retinalExtentDeg = retinalExtentMm/mmPerDeg;
fprintf('Retinal extent mm: %0.1f, degrees (approx): %0.1f\n',retinalExtentMm,retinalExtentDeg);

%% Find spread (analogous to PSF) of each ray bundle where it hits the ey
for hh = 1:length(rayHorizStartPosMm)
    retinalPSFHomologue(hh) = nanstd(retinalHorizPositions(hh,:));
    numberOfRetinalRays(hh) = length(find(~isnan(retinalHorizPositions(hh,:))));
    fprintf('Bundle %d spread (std): %0.3f mm (%0.3f arcmin), based on %d rays\n',hh,retinalPSFHomologue(hh),60*retinalPSFHomologue(hh)/mmPerDeg,numberOfRetinalRays(hh))
end

%% Spatial sampling calcs
pixelsPerDeg = totalPixels/retinalExtentDeg;

%% Max spatial frequcency sampling
maxSpatialFrequencyCpd = 30;
maxSpatialFrequencyDegPerCycle = 1/maxSpatialFrequencyCpd;
maxSpatialFrequencyPixelsPerCycle = pixelsPerDeg*maxSpatialFrequencyDegPerCycle;
fprintf('At %0.1f cycles per deg, have %0.3f deg per cycle, %0.1f pixels per cycle\n', ...
    maxSpatialFrequencyCpd,maxSpatialFrequencyDegPerCycle,maxSpatialFrequencyPixelsPerCycle);

%% Min spatial frequency sampling
minSpatialFrequencyCpd = 1;
minSpatialFrequencyDegPerCycle = 1/minSpatialFrequencyCpd;
minSpatialFrequencyPixelsPerCycle = pixelsPerDeg*minSpatialFrequencyDegPerCycle;
fprintf('At %0.1f cycles per deg, have %0.1f cycles per image\n', ...
    minSpatialFrequencyCpd,retinalExtentDeg/minSpatialFrequencyDegPerCycle);

%% Draw etc.
figure(1);
title(sprintf('Focal length %0.1f, lens diameter %0.1f mm, artificial pupil %0.1f mm', ...
    lensFocalLengthMm,lensDiameterMm,artificialPupilRadiusMm));
drawnow;


