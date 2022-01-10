% Subclass of @Calibrator based on PsychImaging-controlled (Psychtoolbox-3) graphics.
%
% 8/05/2014  npc   Wrote it.
%

classdef SACCPsychImagingCalibrator < Calibrator  
    % Public properties (specific to the @PsychImagingCalibrator class) 
    properties

    end

    % --- PRIVATE PROPERTIES ----------------------------------------------
    properties (Access = private) 
        % handle to screen to be calibrated
        masterWindowPtr;
        
        % handle to the other screen (if it exists)
        slaveWindowPtr;
        
        % array with all the open textures
        texturePointers = [];
        
        % screenRect of screen to be calibrated
        screenRect;
        
        % the original LUT (to be restored upon termination)
        origLUT;
        
        % normalMode. True if LEDs should be in normal mode.  False
        % otherwise.
        normalMode = true;
        
        % logical to physical mapping
        logicalToPhysical = [0:15];
        
        % number of subprimaries
        nSubprimaries = 16;
        
        % number of projector primaries
        nPrimaries = 3;
    end
    
    
    % Public methods
    methods
        % Constructor
        function obj = SACCPsychImagingCalibrator(varargin)  
            % Call the super-class constructor.
            obj = obj@Calibrator(varargin{:});
            
            obj.graphicsEngine = 'PsychImaging';
            
            % Verify validity of screen params values
            obj.verifyScreenParamValues();
        end
    end % Public methods

    % Implementations of required -- Public -- Abstract methods defined in the @Calibrator interface   
    methods
        % Method to set the initial state of the displays
        setDisplaysInitialState(obj, userPrompt);

        % Method to update the stimulus and conduct a single radiometric measurement by 
        % calling the corresponding method of the attached @Radiometer object.
        [measurement, S] = updateStimulusAndMeasure(obj, bgSettings, targetSettings, useBitsPP);

        % Method to ensure that the parameters of the screen match those specified by the user
        obj = verifyScreenParamValues(obj);
        
        % Method to shutdown the Calibrator
        obj = shutdown(obj);    
    end % Implementations of required -- Public -- Abstract methods defined in the @Calibrator interface

    % Private methods that only the PsychImagingCalibrator object can call
    methods (Access = private)  
        
        % Method to change the background and target color
        updateBackgroundAndTarget(obj, bgSettings, targetSettings, useBitsPP);   
    end  % Private methods 

end