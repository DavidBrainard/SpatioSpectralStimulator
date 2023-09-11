function [gaborImageObject] = MakeImageSettingsFromPtCld(ptCldObject,screenCalObj,standardGaborCalObject,screenBgExcitations,stimulusN,options)
% Make a gabor image from the point cloud object.
%
% Syntax:
%    [gaborImageObject] = MakeImageSettingsFromPtCld(ptCldObject,screenCalObj,standardGaborCalObject,screenBgExcitations,stimulusN)
%
% Description:
%    This makes a gabor image with a desired contrast using the point cloud
%    method.
%
% Inputs:
%    ptCldObject               - Structure with the contrasts for all
%                                possible settings using the point cloud
%                                method.
%    screenCalObj              - Screen calibration object.
%    standardGaborCalObject    - Structure with the gabor contrasts
%                                and settings in a cal format.
%    screenBgExcitations       - Screen background cone excitations.
%    stimulusN                 - The size of the stimulus (gabor image) in
%                                pixels.
%
% Outputs:
%    gaborImageObject          - Structure with the gabor contrast image in
%                                image format.
%
% Optional key/value pairs:
%    lightVer                  - Deafult to true. Print out less variables
%                                saved in the structure. It does not affect
%                                making final gabor images, but saving some
%                                time and memory.
%    verbose                   - Boolean. Default true. Controls
%                                plotting and printout.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  dhb,gka,smo      - Wrote it.
%   01/24/22  smo              - Made it work.
%   01/31/22  smo              - It is possible to work on multiple
%                                target contrast gabors inside this
%                                function.
%   02/08/22  smo              - Added an option to print out less variable
%                                saved in the final structure.
%   05/09/22  smo              - Added an option to make a phase shift on
%                                sine image.

%% Set parameters.
arguments
    ptCldObject
    screenCalObj
    standardGaborCalObject
    screenBgExcitations
    stimulusN
    options.lightVer (1,1) = true
    options.verbose (1,1) = true
end

%% Get image from point cloud in cal format.
%
% We want this routine to take contrast explicitly, expressed relative to
% max contrast we set up, when it makes the image.  We will call this
% multiple times to make stimuli of different contrasts.
nContrastPoints = size(standardGaborCalObject.desiredContrastGaborCal,2);
nPhaseShifts = size(standardGaborCalObject.desiredContrastGaborCal,1);

for ss = 1:nPhaseShifts
    for cc = 1:nContrastPoints
        uniqueQuantizedSettingsGaborCal = SettingsFromPointCloud(ptCldObject.contrastPtCld,...
            standardGaborCalObject.desiredContrastGaborCal{ss,cc},ptCldObject.ptCldSettingsCal);
        
        % Print out min/max of settings
        if (options.verbose)
            fprintf('Gabor image min/max settings: %0.3f, %0.3f\n',min(uniqueQuantizedSettingsGaborCal(:)), max(uniqueQuantizedSettingsGaborCal(:)));
        end
        
        % Get contrasts we think we have obtianed
        uniqueQuantizedExcitationsGaborCal = SettingsToSensor(screenCalObj,uniqueQuantizedSettingsGaborCal);
        uniqueQuantizedContrastGaborCal = ExcitationsToContrast(uniqueQuantizedExcitationsGaborCal,screenBgExcitations);
        
        % Plot of how well point cloud method does in obtaining desired contrats
        if (options.verbose)
            figure; clf;
            desiredContrastGaborCal = cell2mat(standardGaborCalObject.desiredContrastGaborCal);
            plot(desiredContrastGaborCal(:),uniqueQuantizedContrastGaborCal(:),'r+');
            axis('square');
            xlabel('Desired L, M or S contrast');
            ylabel('Predicted L, M, or S contrast');
            title('Quantized unique point cloud image method');
        end
        
        % Convert representations we want to take forward to image format. Also, save the results in a structure.
        if(~options.lightVer)
            gaborImageObject.uniqueQuantizedContrastGaborImage{ss,cc} = CalFormatToImage(uniqueQuantizedContrastGaborCal,stimulusN,stimulusN);
            gaborImageObject.uniqueQuantizedSettingsGaborImage{ss,cc} = CalFormatToImage(uniqueQuantizedSettingsGaborCal,stimulusN,stimulusN);
            gaborImageObject.desiredContrastGaborImage{ss,cc} = CalFormatToImage(standardGaborCalObject.desiredContrastGaborCal{ss,cc},stimulusN,stimulusN);
            gaborImageObject.standardPredictedContrastImage{ss,cc} = CalFormatToImage(standardGaborCalObject.standardPredictedContrastGaborCal{ss,cc},stimulusN,stimulusN);
        end
        gaborImageObject.standardSettingsGaborImage{ss,cc} = CalFormatToImage(standardGaborCalObject.standardSettingsGaborCal{ss,cc},stimulusN,stimulusN);
    end
end
end