function [ISETBioDisplayObject,screenCalObjFromISETBio] = MakeISETBioDisplayObj(screenCalObj,screenDistanceMeters,...
    screenSizeMeters,screenSizePixels,T_cones,S,screenGammaMethod,options)
% Make ISETBio display object from the calibration object.
%
% Syntax:
%    [ISETBioDisplayObject] = MakeISETBioDisplayObj(screenCalObj,screenDistanceMeters,...
%     screenSizeMeters,screenSizePixels,T_cones,S,screenGammaMethod)
%
% Description:
%    This makes the ISETBio display objects from the calibration objects
%    with the information of its viewing conditions.
%
%    Note that you don't need the variables of 'T_cones', 'S', and
%    'screenGammaMethod' for making ISETBio object, but we want to make
%    sure if we can get the same values for the reverse calculations (from
%    ISETBio to calibration object). Thus, we added those to make sure if
%    every calculation is correct.
%
% Inputs:
%    screenCalObj            - Target screen calibration object.
%    screenDistanceMeters    - Screen distance from the observer in meters.
%    screenSizeMeters        - Screen size in meters. It should be in the
%                              format [horizontal vertical].
%    screenSizePixels        - Screen resolution in pixels. It should be in
%                              the format [horizontal vertical].
%    T_cones                 - Spectral cone sensitivity function. This is
%                              not needed for making ISETBio object, but to
%                              check if everything is correct when we
%                              calculate backwards from ISETBio to
%                              calibration object.
%    S                       - Spectral wavelength range in cal format.
%    screenGammaMethod       - Screen gamma fitting method.
%
% Outputs:
%    ISETBioDisplayObject    - Created ISETBio display object.
%    screenCalObjFromISETBio - Screen calibration object calculated from
%                              the ISETBio object. This should be the same
%                              with the input screen calibration object.
%
% Optional key/value pairs:
%    verbose                 - Default true. Print out status message if it
%                              is set to true.
%
% See also:
%    SpectralCalCompute, SpectralCalCheck, SpectralCalAnalyze,
%    SpectralCalISETBio

% History:
%   01/21/22  smo            - Wrote it

%% Set parameters.
arguments
    screenCalObj
    screenDistanceMeters (1,1)
    screenSizeMeters (1,2)
    screenSizePixels (1,2)
    T_cones
    S
    screenGammaMethod (1,1)
    options.verbose (1,1) = true
end

%% Create ISETBio display.
%
% Set extra cal data.
extraCalData = ptb.ExtraCalData;
extraCalData.distance = screenDistanceMeters;

% Set screen cal object.
metersToMillimeters = 1000;
screenCalStruct = screenCalObj.cal;
screenCalStruct.describe.displayDescription.screenSizeMM = metersToMillimeters * screenSizeMeters;
screenCalStruct.describe.displayDescription.screenSizePixel = screenSizePixels;
screenCalName = screenCalObj.cal.describe.monitor;

% Make a ISETBio display object here.
ISETBioDisplayObject = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct(screenCalName, screenCalStruct, extraCalData, false);
ISETBioDisplayObject = rmfield(ISETBioDisplayObject,'dixel');

%% Get calibration structure back out to check if every value is correct.
%
% This should match screenCalObj, and we should be able to get same image
% data from the ISETBio scene as from the PTB routines.
screenCalStructFromISETBio = ptb.GeneratePTCalStructFromIsetbioDisplayObject(ISETBioDisplayObject);
screenCalObjFromISETBio = ObjectToHandleCalOrCalStruct(screenCalStructFromISETBio);
SetSensorColorSpace(screenCalObjFromISETBio, T_cones, S);
SetGammaMethod(screenCalObjFromISETBio, screenGammaMethod);

% Let's check that what comes back is what went in.
calCheckList = {'S', 'P_device', 'gammaInput', 'gammaTable', 'P_ambient'};
nCalCheckList = size(calCheckList,2);
for cc = 1:nCalCheckList
    calCheckStr = calCheckList{cc};
    if (any(any(screenCalObj.get(calCheckStr) ~= screenCalObjFromISETBio.get(calCheckStr))))
        error(append(calCheckStr,' is distorted in and out of ISETBio display object'));
    end
end

%% Print out if everything went well.
if (options.verbose)
    disp('ISETBio display object has been successfully created!');
end

end