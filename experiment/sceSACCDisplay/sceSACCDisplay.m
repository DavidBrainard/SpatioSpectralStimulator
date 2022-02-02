function dataOut = sceSACCDisplay(sceneEngineOBJ,testContrast,sceneParamsStruct)
% Compute function for generating a sequence of scenes depicting a
% temporal modulation of a uniform field.
%
% Syntax:
%   dataOut = sceSACCDisplay(obj,testContrast,sceneParamsStruct);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object.  This compute function is set up to return the information
%    required to run a psychophysical experiment on our SACC display, and
%    also so that can be used within the ISETBioCSFGenerator framework.
%
%    Because the computations to generate scenes for each contrast are
%    relatively slow, we predefine the allowable contrasts.  This is OK,
%    because we can limit the contrasts that we'll call for to the same
%    ones for use by Quest+ in the CSF generator code.
%
%    Thus in the implementation, we actually just pass in as part of the
%    sceneParameters structure the predefined scenes that we will use, and
%    just look up the answer and pass it back.
%
%    This routine uses the statusReport field of the returned structure to
%    hold the RGB data that we want.  The status report is an extra field
%    understood by the sceneEngine object, that can be used for returning
%    custom data.  We're overloading that a bit by returning something
%    substantive, but we don't think that does any harm.
%
%    Also note that because of the extensive precomputation required, the
%    default parameters returned by this routine are a skeleton. If you
%    passed the defaults in, you'd get an error because you would not have
%    any precomputed scenes or RGB images.  We think this is also OK,
%    because we can simply avoid ever doing this in our calling code for
%    this particular scene engine.
% 
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%               dataOut = sceUniformFieldTemporalModulation();
%           it does not compute anything and simply returns a struct with the 
%           defaultParams that define the scene.
%
%       [2] If called with arguments, as it is from a parent @sceneEngine object,
%           it computes a cell array of scenes defining the frames of a
%           stimulus and the temporal support of the frames. These are
%           returned as named fields of the returned dataOut struct.
%
%    All scene functions used with the sceneEngine class must conform to
%    this API.
%
% Inputs:
%    sceneEngineOBJ             - Calling @sceneEngine object.  This is
%                                  not currently used, but passing it allows us
%                                  flexibility in the future and matches
%                                  conventions for the other classes in
%                                  this toolbox.
%    testContrast                - Scalar providing the contrast for the
%                                  scene to be generated.                           
%    sceneParamsStruct           - Struct containing properties of the
%                                  scene understood by this function.
%                                  As noted above, execute
%                                  sceSACCDisplay at the
%                                  command line to see the structure's
%                                  fields and default values.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams that define the scene
%
%             - If called from a parent @sceneEngine, the returned
%               struct is organized as follows:
%                 .sceneSequence : a cell array of scenes defining the frames of the generated grating scene sequence                            
%                 .temporalSupport : the temporal support of the frames of the generated grating scene sequence, in seconds
%                 .statusReport : a structure containing additional
%                                 information, and particularly the RGB images we need to show the stimulus
%                                 at each predefined contrast.
%
% Optional key/value input arguments:
%    None.
%
% Examples:
%    The source code contains examples.
%
% See Also:
%     t_sceneGeneration, t_thresholdEngine

% History:
%    01/25/22  dhb, smo  Started on this.
%    01/26/22  dhb, smo  Wrote version 1.

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Check that sceneParameters contains at least one precomputed contrast
    if (isempty(sceneParamsStruct.predefinedContrasts))
        error('Passed sceneParameters struct contains no predefined contrast data');
    end

    % Called for real, so we need to return the information for the passed
    % contrast.
    whichContrast = find(testContrast == sceneParamsStruct.predefinedContrasts);
    if (isempty(whichContrast))
        error('Contrast requested is not in predefined list');
    end

    % Pick out precomputed items for return
    dataOut.sceneSequence = sceneParamsStruct.predefinedSceneSequences{whichContrast};
    dataOut.temporalSupport = sceneParamsStruct.predefinedTemporalSupport;
    dataOut.statusReport.RGBimage = sceneParamsStruct.predefinedRGBImages{whichContrast};
end

function p = generateDefaultParams()
    p = struct(...
        'predefinedContrasts', [], ...              % This can only generate scenes for these contrasts
        'predefinedRGBImages', {}, ...              % Precomputed RGB images for each contrast
        'predefinedSceneSequences', {}, ...         % Precomputed scenes for each contrast           
        'predefinedTemporalSupport', [] ...         % Common precomputed temporal support
    );
end
