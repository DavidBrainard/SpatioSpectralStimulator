%% ISETBioCalcsMel.m
%
% Use ISETBio to compute cone splatter given an ISETBio scene
% that describes the image.
%
% The scene is most naturally set up using SpectralCalISETBioMel.
% That script stores the scene and other information in a .mat
% file following a naming convention it sets up based on the
% parameters, and this routine finds the file using the same
% naming convention.
%
% There are examples in the header comments below that set this
% script loose over various parameters, given that the scenes
% have been set up.  These illustrate how to explore individual
% variation, field location, and effects of LCA as a function
% of spatial frequency.

%% TODO
%
% Consider turning off the off axis attenuation. Shouldn't
% affect contrast much, and complicates correction. 
%
% Compute with noLCA to check that contrasts are sf independent
% in that case.  DONE.
%
% Provide method to deconvolve an image by the polychromatic CSF.
% Veryify that this removes effects of CA on splatter.
%
% Save and allow reuse of OI for de-convolving, so we can ask
% what happens if we correct for one SF and then apply at 
% others.
%
% Understand how well we can approximate a decovolved image
% with the primaries that produced the raw image.
%
% Look at effects of individual variation in LCA, and possibly
% TCA on all this.

% Asano parameter ranges.
%
% Asano et al. give the following population SD's for the individual
% difference parameters (their Table 5, Step 2 numbers:
%   Lens    - 18.7%
%   Macular - 36.5%
%   L Density - 9%
%   M Density - 9%
%   S Density - 7.4%
%   L Shift   - 2 nm
%   M Shift   - 1.5 nm
%   S Shift   - 1.3 nm
% Knowing these can be useful in choosing individual difference
% parameter ranges to explore.

%% Run out sf example
%{
ion

%}
 l
%% Run out MP variation example
%{
clear;
conditionNameList = {'MelDirected1'}; % {'MelDirected1', 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2];
gaborSdDeg = 100;
stimulusSizeDeg = 4;
stimConeEccDeg = 25;

fieldSizeDeg = 5;
eccXDegList = -5; %eccXDegList = [0 -5 -10 -15 -20];
aoRender = false;
noLCA = false;
ageISETBioList = 32; %[20 32 60];
dLensISETBioList = 0; %[-18.7 0 18.7];
dMacISETBioList = [1.5*-36.5 -36.5 -36.5/2 0 36.5/2 36.5 1.5*36.5];

for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        for ee = 1:length(eccXDegList)
            for aa = 1:length(ageISETBioList)
                for ll = 1:length(dLensISETBioList)
                    for mm = 1:length(dMacISETBioList)
                        lmsConeContrast(:,cc,ss,ee,aa,ll,mm) = ISETBioCalcsMel(conditionNameList{cc},stimConeEccDeg,sineFreqCyclesPerDegList(ss), ...
                                gaborSdDeg,stimulusSizeDeg, ...
                                fieldSizeDeg, ...
                                eccXDegList(ee), ...
                                aoRender, noLCA, ...
                                ageISETBioList(aa), dLensISETBioList(ll), dMacISETBioList(mm) ...
                            );
                    end
                end
            end
        end
    end
end

% Plot of results
lmsContrastByDMac = squeeze(lmsConeContrast(:,:,:,:,:,:,1:end));
figure; clf; hold on;
plot(dMacISETBioList,100*lmsContrastByDMac(1,:),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(dMacISETBioList,100*lmsContrastByDMac(2,:),'go','MarkerFaceColor','g','MarkerSize',12);
plot(dMacISETBioList,100*lmsContrastByDMac(3,:),'bo','MarkerFaceColor','b','MarkerSize',12);
plot(dMacISETBioList,100*lmsContrastByDMac(1,:),'r','LineWidth',2);
plot(dMacISETBioList,100*lmsContrastByDMac(2,:),'g','LineWidth',2);
plot(dMacISETBioList,100*lmsContrastByDMac(3,:),'b','LineWidth',2);
set(gca,'FontName','Helvetica','FontSize',18);
xlabel('Macular Pigment Density Adj (Percent)','FontName','Helvetica','FontSize',20);
ylabel('Contrast (percent)','FontName','Helvetica','FontSize',20);
legend({'L cones','M cones','S cones'},'FontName','Helvetica','FontSize',14,'Location','NorthWest');
ylim([0 10]);

%}

%% Run out lens variation example
%{
clear;
conditionNameList = {'MelDirected1'}; % {'MelDirected1', 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2];
gaborSdDeg = 100;
stimulusSizeDeg = 4;
stimConeEccDeg = 25;

fieldSizeDeg = 5;
eccXDegList = -5; %eccXDegList = [0 -5 -10 -15 -20];
aoRender = false;
noLCA = false;
ageISETBioList = 32; %[20 32 60];
dLensISETBioList = [1.5*-18.7 -18.7 -18.7/2 0 18.7/2 18.7 1.5*18.7];
dMacISETBioList = 0;
for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        for ee = 1:length(eccXDegList)
            for aa = 1:length(ageISETBioList)
                for ll = 1:length(dLensISETBioList)
                    for mm = 1:length(dMacISETBioList)
                        lmsConeContrast(:,cc,ss,ee,aa,ll,mm) = ISETBioCalcsMel(conditionNameList{cc},stimConeEccDeg,sineFreqCyclesPerDegList(ss), ...
                                gaborSdDeg,stimulusSizeDeg, ...
                                fieldSizeDeg, ...
                                eccXDegList(ee), ...
                                aoRender, noLCA, ...
                                ageISETBioList(aa), dLensISETBioList(ll), dMacISETBioList(mm) ...
                            );
                    end
                end
            end
        end
    end
end

% Plot of results
lmsContrastBydLensI = squeeze(lmsConeContrast(:,:,:,:,:,1:end,:));
figure; clf; hold on;
plot(dLensISETBioList,100*lmsContrastBydLensI(1,:),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(dLensISETBioList,100*lmsContrastBydLensI(2,:),'go','MarkerFaceColor','g','MarkerSize',12);
plot(dLensISETBioList,100*lmsContrastBydLensI(3,:),'bo','MarkerFaceColor','b','MarkerSize',12);
plot(dLensISETBioList,100*lmsContrastBydLensI(1,:),'r','LineWidth',2);
plot(dLensISETBioList,100*lmsContrastBydLensI(2,:),'g','LineWidth',2);
plot(dLensISETBioList,100*lmsContrastBydLensI(3,:),'b','LineWidth',2);
set(gca,'FontName','Helvetica','FontSize',18);
xlabel('Lens Pigment Density Adj (Percent)','FontName','Helvetica','FontSize',20);
ylabel('Contrast (percent)','FontName','Helvetica','FontSize',20);
legend({'L cones','M cones','S cones'},'FontName','Helvetica','FontSize',14,'Location','NorthWest');
ylim([0 10]);

%}

%% Run out eccentricity example
%{
clear;
conditionNameList = {'MelDirected1'}; % {'MelDirected1', 'IsochromaticControl'};
sineFreqCyclesPerDegList = [0.2];
gaborSdDeg = 100;
stimulusSizeDeg = 4;
stimConeEccDeg = 20;

fieldSizeDeg = 5;
eccXDegList = [0 -2.5 -5 -10 -15 -20];
aoRender = false;
noLCA = false;
ageISETBioList = 32; %[20 32 60];
dLensISETBioList = 0; %[1.5*-18.7 -18.7 -18.7/2 0 18.7/2 18.7 1.5*18.7];
dMacISETBioList = 0;
for cc = 1:length(conditionNameList)
    for ss = 1:length(sineFreqCyclesPerDegList)
        for ee = 1:length(eccXDegList)
            for aa = 1:length(ageISETBioList)
                for ll = 1:length(dLensISETBioList)
                    for mm = 1:length(dMacISETBioList)
                        lmsConeContrast(:,cc,ss,ee,aa,ll,mm) = ISETBioCalcsMel(conditionNameList{cc},stimConeEccDeg,sineFreqCyclesPerDegList(ss), ...
                                gaborSdDeg,stimulusSizeDeg, ...
                                fieldSizeDeg, ...
                                eccXDegList(ee), ...
                                aoRender, noLCA, ...
                                ageISETBioList(aa), dLensISETBioList(ll), dMacISETBioList(mm) ...
                            );
                    end
                end
            end
        end
    end
end

% Plot of results
lmsContrastByEccX = squeeze(lmsConeContrast(:,:,1:end,:,:,:,:));
figure; clf; hold on;
plot(abs(eccXDegList),100*lmsContrastByEccX(1,:),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(abs(eccXDegList),100*lmsContrastByEccX(2,:),'go','MarkerFaceColor','g','MarkerSize',12);
plot(abs(eccXDegList),100*lmsContrastByEccX(3,:),'bo','MarkerFaceColor','b','MarkerSize',12);
plot(abs(eccXDegList),100*lmsContrastByEccX(1,:),'r','LineWidth',2);
plot(abs(eccXDegList),100*lmsContrastByEccX(2,:),'g','LineWidth',2);
plot(abs(eccXDegList),100*lmsContrastByEccX(3,:),'b','LineWidth',2);
set(gca,'FontName','Helvetica','FontSize',18);
xlabel('Eccentricity (deg)','FontName','Helvetica','FontSize',20);
ylabel('Contrast (percent)','FontName','Helvetica','FontSize',20);
legend({'L cones','M cones','S cones'},'FontName','Helvetica','FontSize',14,'Location','NorthWest');
ylim([0 10]);

%}

%% Baseline condition test
%{
clear;
conditionName = 'MelDirected1';
sineFreqCyclesPerDeg = 0.2;
gaborSdDeg = 100;
stimulusSizeDeg = 4;
stimConeEccDeg = 25;

fieldSizeDeg = 5;
eccXDeg = -5;
aoRender = false;
noLCA = false;
ageISETBio = 32;
dLensISETBio = 0;
dMacISETBio = 0;

[lmsConeContrast] = ISETBioCalcsMel(conditionName,stimConeEccDeg,sineFreqCyclesPerDeg, ...
        gaborSdDeg,stimulusSizeDeg, ...
        fieldSizeDeg, ...
        eccXDeg, ...
        aoRender, noLCA, ...
        ageISETBio, dLensISETBio, dMacISETBio ...
    );

%}

% ISETBioCalcsMel
function [lmsConeContrast] = ISETBioCalcsMel(conditionName,stimConeEccDeg,sineFreqCyclesPerDeg, ...
            gaborSdDeg,stimulusSizeDeg, ...
            fieldSizeDeg, ...
            eccXDeg, ...
            aoRender,noLCA, ...
            ageISETBio, dLensISETBio, dMacISETBio ...
        )

%% Initialize
close all;

%% Boring
legendFontSize = 12;
labelFontSize = 18;
titleFontSize = 20;
axisFontSize = 16;

%% Parameter specifications for scene that are not passed.
screenGammaMethod = 2;
sceneInputStr = sprintf('%s_StimConeEcc_%0.1f_Size_%0.1f_Sf_%0.1f_Sd_%0.1f_GammaMethod_%d', ...
    conditionName,stimConeEccDeg,stimulusSizeDeg,sineFreqCyclesPerDeg,gaborSdDeg,screenGammaMethod);

% Load the scene data according to parameters above
projectFiledir = getpref('SpatioSpectralStimulator','SACCMelanopsin');
sceneInputSubdir = fullfile(projectFiledir,sceneInputStr);
sceneInputFilename = fullfile(sceneInputSubdir,'sceneOutput.mat');
theData = load(sceneInputFilename);
disp('Data loaded');

%% Mosaic settings not so likely to vary
eccYDeg = 0;
subjectID = 6;
mosaicOutputStr = sprintf('%s_StimConeEcc_%0.1f_EccX_%0.1f_EccY_%0.1f_FieldSize_%0.1f_AO_%d_NoLCA_%d_Age_%d_macAdj_%0.1f_lensAdj_%0.1f', ...
    conditionName,stimConeEccDeg,eccXDeg,eccYDeg,fieldSizeDeg,aoRender,noLCA,ageISETBio,dMacISETBio,dLensISETBio);

%% Use ISETPipelineToolbox wrapper as a way to the isetbio computations
%
% These parameters define key addtional properties we are likely to vary
% Not likely to vary these in the short run, so not coded in output name for now
randSeed = false;                              % False means, no randomness
eccVars = false;
pupilDiamMM = 3;
defocusDiopters = 0;
zernikeDataBase = 'Polans2015';

%% Setup output names and places
outputDir = fullfile(sceneInputSubdir,mosaicOutputStr);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
mosaicOutputFilename = fullfile(outputDir,[sprintf('mosaicOutput') '.mat']);

%% Compute radial eccentricity
eccDeg = sqrt(eccXDeg^2 + eccYDeg^2);

%% Generate the photopigment, macular pigment density, and lens density to use
%
% What we're going to do is override the oi and cMosaic's view of what
% these are, so we can control them completely.  We'll calculate these
% first so that we can insert the cMosaic parameters on creation, which is
% easier than doing it ex post.  See tutorial t_conesPhotopigment.
%
% Lens to use
lensTransmittance = LensTransmittance(theData.colorDirectionParams.wls,'Human','CIE', ageISETBio, pupilDiamMM);
unadjustedLensDensity = -log10(lensTransmittance);
adjustedLensDensity = unadjustedLensDensity * (1 + dLensISETBio/100);
lensTransmittance = 10.^-adjustedLensDensity;
lensObject = Lens('wave',theData.colorDirectionParams.wls,'unitDensity',-log10(lensTransmittance),'density',1);

% Plot of lens density versus age
plotAges = linspace(20,80,100);
for aa = 1:length(plotAges)
    lensTransmittanceAA = LensTransmittance(theData.colorDirectionParams.wls,'Human','CIE', plotAges(aa), pupilDiamMM);
    unadjustedLensDensityAA = -log10(lensTransmittanceAA);
    plotLensDensity(aa) = max(unadjustedLensDensityAA(:));
end
figure; clf; hold on
plot(plotAges,plotLensDensity,'b','LineWidth',6);
xlabel('Age (years)','FontName','Helvetica','FontSize',labelFontSize);
ylabel('Lens Peak Density','FontName','Helvetica','FontSize',labelFontSize);
ylim([0 6]);
saveas(gcf,fullfile(outputDir,'lensDensityWithAge.tiff'),'tiff');

%% Macular transmittance for the specified field size
macTransmittance = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',2*eccDeg);
unadjustedMacDensity = -log10(macTransmittance);
adjustedMacDensity = unadjustedMacDensity * (1 + dMacISETBio/100);
macTransmittance = 10.^-adjustedMacDensity;
macObject = Macular('wave',theData.colorDirectionParams.wls,'unitDensity',-log10(macTransmittance),'density',1);

% Nominal macular tranmittance for explanatory plot
macTransmittanceNominal = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',2*eccDeg);
macTranFig = figure; clf; hold on;
plot(theData.colorDirectionParams.wls,macTransmittanceNominal,'b','LineWidth',6);
xlabel('Wavelength (nm)','FontName','Helvetica','FontSize',labelFontSize);
ylabel('Macular Pigment Transmittance','FontName','Helvetica','FontSize',labelFontSize);
ylim([0 1]);
saveas(gcf,fullfile(outputDir,'macTransmittance.tiff'),'tiff');

% Plot of macular pigment density versus eccentricity
plotEccs = linspace(0,25,100);
for ee = 1:length(plotEccs)
    macTransmittanceEE = MacularTransmittance(theData.colorDirectionParams.wls,'Human','CIE',2*plotEccs(ee));
    unadjustedMacDensityEE = -log10(macTransmittanceEE);
    plotMacDensity(ee) = max(unadjustedMacDensityEE(:));
end
figure; clf; hold on
plot(plotEccs,plotMacDensity,'b','LineWidth',6);
xlabel('Eccentricity (degs)','FontName','Helvetica','FontSize',labelFontSize);
ylabel('Macular Pigment Peak Density','FontName','Helvetica','FontSize',labelFontSize);
ylim([0 max(plotMacDensity(:))]);
saveas(gcf,fullfile(outputDir,'macDensityWithEccentricty.tiff'),'tiff');

%% Photopigment absorptance
%
% The individual differences model allows for shifting the absorbance along
% the wavelength axis. Do that here.
asanoConeParams = DefaultConeParams('cie_asano');
asanoConeParams.fieldSizeDegrees = 2*eccDeg;
asanoConeParams.ageYears = ageISETBio;
asanoConeParams.pupilDiamMM = pupilDiamMM;
asanoConeParams.indDiffParams = theData.colorDirectionParams.psiParamsStruct.coneParams.indDiffParams;
[~,~,~,~,cieConeParams,cieStaticParams] = ...
    ComputeCIEConeFundamentals(MakeItS(theData.colorDirectionParams.wls),asanoConeParams.fieldSizeDegrees,...
    asanoConeParams.ageYears,asanoConeParams.pupilDiamMM,[],[],[], ...
    [],[],[],asanoConeParams.indDiffParams);

% Photopigment absorbance
temp = load('T_log10coneabsorbance_ss');
photopigmentAbsorbance = 10.^SplineCmf(temp.S_log10coneabsorbance_ss,temp.T_log10coneabsorbance_ss,theData.colorDirectionParams.wls,2);
clear temp

% Get axial density.  This scales the normalized photopigment absorbance we
% computed just above.
axialDensity = cieConeParams.axialDensity .* (1 + asanoConeParams.indDiffParams.dphotopigment/100);

% Get OS length so we can compute axial density from first principles
%
% Rodieck, The First Steps in Seeing, p. 472 gives a value for
% axial specific density for rods and cones as about 0.015 /um.
[osLengthMicrons, osLengthMicronsFoveal] = cMosaic.outerSegmentLengthFromEccentricity(eccDeg);
specificDensity = 0.015;
axialDensityFromOSLength = specificDensity*osLengthMicrons;
fprintf('Axial density at %d degs from CIE: %0.3f, %0.3f, %0.3f; from OS length and specific density %0.3f\n', eccDeg, ...
    cieConeParams.axialDensity(1),cieConeParams.axialDensity(2),cieConeParams.axialDensity(3),axialDensityFromOSLength);

% Set up photopigment object
%
% This is based on the axialDensity (parameter name 'opticalDensity'),
% absorbance (parameter name 'absorbance'), and quantalEfficiency
% (parameter name 'peakEfficiency'). Someday we will rewrite so that naming 
% conventions are consistent across different places we implement them, perhaps.
% 
% ISETBio wants absorbance in columns, while PTB wants it in rows, so we
% also need to transposes.
photopigmentObject = cPhotoPigment('wave', theData.colorDirectionParams.wls,...
    'opticalDensity',axialDensity,'absorbance',photopigmentAbsorbance', ...
    'peakEfficiency',cieStaticParams.quantalEfficiency );

% Get the cone fundamentals we think we are using
T_cones_ISETBio = photopigmentObject.energyFundamentals';
for cc = 1:3
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) .* lensTransmittance;
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) .* macTransmittance;
    T_cones_ISETBio(cc,:) = T_cones_ISETBio(cc,:) / max(T_cones_ISETBio(cc,:));
end
fundamentalsFig = figure; clf; hold on;
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(1,:),'r','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(2,:),'g','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(3,:),'b','LineWidth',6);
xlabel('Wavelength (nm)');
ylabel('Exication Probability');
figure(fundamentalsFig);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(1,:),'k-','LineWidth',3);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(2,:),'k-','LineWidth',3);
plot(theData.colorDirectionParams.wls,T_cones_ISETBio(3,:),'k-','LineWidth',3);

fundamentalsWithMelBaselineFig = figure; clf; hold on;
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(1,:),'r','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(2,:),'g','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_cones(3,:),'b','LineWidth',6);
plot(theData.colorDirectionParams.wls,theData.colorDirectionParams.T_receptors(4,:),'c','LineWidth',6);
legend({'L', 'M', 'S', 'Mel'},'FontName','Helvetica','FontSize',12);
xlabel('Wavelength (nm)','FontName','Helvetica','FontSize',labelFontSize);
ylabel('Exication Probability','FontName','Helvetica','FontSize',labelFontSize);
saveas(fundamentalsWithMelBaselineFig ,fullfile(outputDir,'fundamentalsWithMelBaselineFig.tiff'),'tiff');

%% Create and setup cone mosaic
%
% Use the ISETPipelineToolbox machinery.

%% AO version allows us to get rid of optics in effect
if (aoRender)
    if (eccVars)
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'subjectID', 0, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', 0, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end

    % Normal optics structure. Allow specified defocus.
else
    if (eccVars)
        % Build normal optics structure.
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDeg, eccYDeg, ...
            'fovealDegree', fieldSizeDeg, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', theData.colorDirectionParams.wls, ...
            'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
            'macular', macObject, ...       % custom macular pigment object
            'pigment', photopigmentObject, ... % custom photopigment object
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end
end

%% Compute optical image
ss = 1; cc = 1;
theConeMosaic.PSF = oiSet(theConeMosaic.PSF,'lens',lensObject);
theOI = oiCompute(theData.ISETBioGaborObject.ISETBioGaborScene{ss,cc}, theConeMosaic.PSF);

% Some plots
sampledWavelengths = [420:40:700];
psfRangeArcMin = 10;
psfFig = visualizePSFsubStack(theOI, sampledWavelengths, psfRangeArcMin);
saveas(psfFig,fullfile(outputDir,'psfFigure.tiff'),'tiff');

uData = oiPlot(theOI, 'irradiance image with grid', [], 40);
saveas(gcf,fullfile(outputDir,'oiIrradiance.tiff'),'tiff');

%% Compute cone responses
theConeResponses = theConeMosaic.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered', ...
    'lowOpticalImageResolutionWarning',true);
%theConeMosaic.Mosaic.plot('excitations',theConeResponses);

%% Visualize the mosaic itself
mosaicFigStruct = theConeMosaic.Mosaic.visualize;
title('Cone Mosaic');
saveas(mosaicFigStruct.figureHandle,fullfile(outputDir,'mosaicImage.tiff'),'tiff');

%% Visualize noise-free responses
%
% Basic activation map
mosaicActivationFigStruct = theConeMosaic.Mosaic.visualize( ...
    'activation', theConeResponses ...
    );
title('Activation Map');
saveas(mosaicActivationFigStruct.figureHandle,fullfile(outputDir,'activationMap.tiff'),'tiff');

%% Compute a contrast responses array
theConeContrasts = zeros(size(theConeResponses));
lConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.lConeIndices));
meanLConeResponse = mean(lConeResponses(:));
theConeContrasts(1,1,theConeMosaic.Mosaic.lConeIndices) = ...
    (theConeResponses(1,1,theConeMosaic.Mosaic.lConeIndices) - meanLConeResponse)/meanLConeResponse;
mConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.mConeIndices));
meanMConeResponse = mean(mConeResponses(:));
theConeContrasts(1,1,theConeMosaic.Mosaic.mConeIndices) = ...
    (theConeResponses(1,1,theConeMosaic.Mosaic.mConeIndices) - meanMConeResponse)/meanMConeResponse;
sConeResponses = squeeze(theConeResponses(1,1,theConeMosaic.Mosaic.sConeIndices));
meanSConeResponse = mean(sConeResponses(:));
theConeContrasts(1,1,theConeMosaic.Mosaic.sConeIndices) = ...
    (theConeResponses(1,1,theConeMosaic.Mosaic.sConeIndices) - meanSConeResponse)/meanSConeResponse;

% Visualize a cut
horizontalSliceYcoordDegs = 0;
horizontalSliceThicknessDegs = 0.08;
coneTypesToVisualize = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
maxResponse = 1.2*max(theConeResponses(:));
% visualizedResponseScalingDegs = 0.5;
% hFig = theConeMosaic.Mosaic.visualizeHorizontalConeActivationProfiles(...
%     theConeResponses, coneTypesToVisualize, ...
%     horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
%     maxResponse, visualizedResponseScalingDegs);
mosaicSliceFig = theConeMosaic.Mosaic.visualizeHorizontalConeActivationProfiles(...
    theConeResponses, coneTypesToVisualize, ...
    horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
    maxResponse,  ...
    'onTopOfConeMosaicActivation', false, ...
    'figureHandle', [], ...
    'axesHandle', [], ...
    'fontSize', 20);
set(gca,'FontName','Helvetica','FontSize',34);
xlabel('Position (degs)','FontName','Helvetica','FontSize',36);
ylabel('Excitation Rate','FontName','Helvetica','FontSize',36);
xlim([eccXDeg-fieldSizeDeg*0.75/2 eccXDeg+fieldSizeDeg*0.75/2]);
saveas(mosaicSliceFig,fullfile(outputDir,'mosaicSliceFig.tiff'),'tiff');

mosaicContrastSliceFig = theConeMosaic.Mosaic.visualizeHorizontalConeActivationProfiles(...
    100*theConeContrasts, coneTypesToVisualize, ...
    horizontalSliceYcoordDegs, horizontalSliceThicknessDegs, ...
    100,  ...
    'onTopOfConeMosaicActivation', false, ...
    'figureHandle', [], ...
    'axesHandle', [], ...
    'fontSize', 20);
set(gca,'FontName','Helvetica','FontSize',34)
xlabel('Position (degs)','FontName','Helvetica','FontSize',36);
ylabel('Contrast (percent)','FontName','Helvetica','FontSize',36);
xlim([eccXDeg-fieldSizeDeg*0.7/2 eccXDeg+fieldSizeDeg*0.7/2]);
ylim([-10 10]);
saveas(mosaicContrastSliceFig,fullfile(outputDir,'mosaicContrastSliceFig.tiff'),'tiff');

%% Choose a circular ROI that captures most of the field and compute contrast
%
% We compute contrast across all cones of each type in the ROI.  Using the
% circular ROI eliminates edge effects.
roiCircle = regionOfInterest('shape','ellipse',...
    'center',[eccXDeg eccYDeg],...
    'majorAxisDiameter',fieldSizeDeg*0.7,...
    'minorAxisDiameter',fieldSizeDeg*0.7);

% Visualize the ROI to make sure it is in a sensible place
theConeMosaic.Mosaic.plot('roi',theConeResponses, 'roi',roiCircle);

% Get L, M, and S cone responses out of ROI.  Assumes just one
% instance/frame, which is what we are running.
lConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','L'));
mConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','M'));
sConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
    'all excitations',theConeResponses, ...
    'cone type','S'));
lConeContrast = (max(lConeResponses(:))-min(lConeResponses(:)))/(max(lConeResponses(:))+min(lConeResponses(:)));
mConeContrast = (max(mConeResponses(:))-min(mConeResponses(:)))/(max(mConeResponses(:))+min(mConeResponses(:)));
sConeContrast = (max(sConeResponses(:))-min(sConeResponses(:)))/(max(sConeResponses(:))+min(sConeResponses(:)));

% Compute residual contrast seen by each class of cones
fprintf('L cone contrast: %0.1f%%; M cone contrast: %0.1f%%; S cone contrast: %0.1f%%\n',100*lConeContrast,100*mConeContrast,100*sConeContrast);
lmsConeContrast = [lConeContrast ; mConeContrast; sConeContrast];

%% Verify that the effective spectral sensitivities in ISETBio match intended
%
% Produce a set of monochromatic images at each sample wavelength.  The photon
% level is scene radiance in photons/sr-m2-nm-sec.
%
% This is slow but checked out well for cases I tried.
VERIFY_ISETBIO_FUNDAMENTALS = false;
if (VERIFY_ISETBIO_FUNDAMENTALS)
    pixelSize = sceneGet(theData.ISETBioGaborObject.ISETBioGaborScene{ss,cc},'rows');
    photonsPerSrM2NMSec = 1e25;
    for ww = 1:length(theData.colorDirectionParams.wls)
        % Set up a dummy scene. Spatially uniform with a black body
        % spectrum of 5000 degK.  We'll replace the scene contents just below.
        scene{ww} = sceneCreate('uniform bb',pixelSize,5000,theData.colorDirectionParams.wls);

        % Use small field of view to minimize effects of eccentricity, and also
        % so we don't need too many pixels (for efficiency in this demo).
        scene{ww} = sceneSet(scene{ww},'fov',stimulusSizeDeg);

        % Get photons and rewrite to be monochromatic constant power in
        % photons/sec-nm.
        photons = sceneGet(scene{ww},'photons');
        photons = zeros(size(photons));
        photons(:,:,ww) = photonsPerSrM2NMSec*ones(pixelSize,pixelSize);
        scene{ww} = sceneSet(scene{ww},'photons',photons);
    end

    %% Compute the retinal image and cone excitations for each scene
    %
    % Use this to get the cone fundamental that ISETBio is effectively
    % using.
    T_quantalFundamentalsFromCMosaic = zeros(3,length(theData.colorDirectionParams.wls));
    T_energyFundamentalsIFromCMosaic = zeros(3,length(theData.colorDirectionParams.wls));

    for ww = 1:length(theData.colorDirectionParams.wls)
        % Compute retinal image
        oiComputed{ww} = oiCompute(scene{ww}, theConeMosaic.PSF);

        % Compute noise free cone excitations
        coneExcitations{ww} = theConeMosaic.Mosaic.compute(oiComputed{ww},'nTrials', 1);

        % Find L, M and S cone excitations at this wavelength.  This is
        % accomplished by extracting the mean response of each cone type from
        % the excitations computed just above and diviting by the input power
        % in the scene.
        %
        % We expect the answer to be proportional to the fundamentals we
        % computed above, because we have not (yet) accounted for the geometry
        % between radiance in the scene and retinal irradiance, nor the cone
        % integration time.
        % Get L, M, and S cone responses out of ROI.  Assumes just one
        % instance/frame, which is what we are running.
        lConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','L'));
        mConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','M'));
        sConeResponses = squeeze(theConeMosaic.Mosaic.excitations('roi',roiCircle,...
            'all excitations',coneExcitations{ww}, ...
            'cone type','S'));
        T_quantalFundamentalsISETBio(1,ww) = mean(lConeResponses(:))/photonsPerSrM2NMSec;
        T_quantalFundamentalsISETBio(2,ww) = mean(mConeResponses(:))/photonsPerSrM2NMSec;
        T_quantalFundamentalsISETBio(3,ww) = mean(sConeResponses(:))/photonsPerSrM2NMSec;
    end

    for ii = 1:3
        T_energyFundamentalsIFromCMosaic(ii,:) = EnergyToQuanta(theData.colorDirectionParams.wls,T_quantalFundamentalsISETBio(ii,:)')';
        T_energyFundamentalsIFromCMosaic(ii,:) = T_energyFundamentalsIFromCMosaic(ii,:)/max(T_energyFundamentalsIFromCMosaic(ii,:));
    end

    figure(fundamentalsFig);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(1,:),'y:','LineWidth',3);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(2,:),'y:','LineWidth',3);
    plot(theData.colorDirectionParams.wls,T_energyFundamentalsIFromCMosaic(3,:),'y:','LineWidth',3);
end
figure(fundamentalsFig);
legend({sprintf('%0.1f%%',100*lConeContrast),sprintf('%0.1f%%',100*mConeContrast),sprintf('%0.1f%%',100*sConeContrast)}, ...
    'FontName','Helvetica','FontSize',12);
saveas(fundamentalsFig,fullfile(outputDir,'fundamentalsFig.tiff'),'tiff');

%% Save out what we need 
colorDirectionParams = theData.colorDirectionParams;
save(mosaicOutputFilename, ...
    'lConeResponses','mConeResponses','sConeResponses', ...
    'lConeContrast','mConeContrast','sConeContrast', 'lmsConeContrast', ...
    'T_cones_ISETBio', 'colorDirectionParams', 'eccXDeg', 'eccYDeg', 'fieldSizeDeg', ...
    'aoRender', 'noLCA', 'ageISETBio', 'dLensISETBio', 'dMacISETBio', ...
    '-v7.3');
disp('Data has been saved successfully!');
