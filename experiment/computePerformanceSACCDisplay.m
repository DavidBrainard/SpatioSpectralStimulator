function [correct] = computePerformanceSACCDisplay(...
            nullRGBImage, testRGBIimage, ...
            theSceneTemporalSupportSeconds,displayControlStruct)
% Run one trial of a psychophysical experiment
%
% Syntax:
%    [correct] = computePerformanceSACCDisplay(...
%            nullRGBImage, testRGBIimage, ...
 %           theSceneTemporalSupportSeconds,displayControlStruct)
%
% Description:
%     Run one trial of a psychophysical experiment and return correct or
%     incorrect.  The trial is TAFC.  The two stimuli are nullRGBImage and
%     testRGBImage.  Subject is correct if he/she chooses testRGBImage.
%
% Inputs:
%     nullRGBImage             - 
%     testRGBImage             - 
%     temporalSupport          - Temporal support vector (in seconds) for
%                                scene sequences.
%     displayControlStruct     - 
%
%
% Outputs:
%     correct                - 1 ifcorrect and 0 if incorrect.
%
% Optional key/value pairs:
%     None.
%
% See also
%   t_thresholdEngine, t_spatialCsf, computeThresholdTAFC
%

% History:
%   10/23/20  dhb  Comments.

% Empty responses
responses = [];

% Train the classifier.
%
% If trainNoiseFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise (typically 'none' or
% 'random') should be used in the training.

if (~isempty(trainNoiseFlag))
    % Generate stimulus for training, NULL stimulus
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        nullScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainNoiseFlag});
    
    % Generate stimulus for training, TEST stimulus
    [inSampleTestStimResponses, responseTemporalSupportSeconds] = theNeuralEngine.compute(...
        testScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainNoiseFlag});
    
    if (visualizeAllComponents)
        if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
            diffResponse = inSampleTestStimResponses(trainNoiseFlag) - inSampleNullStimResponses(trainNoiseFlag);
            % Visualize the activation
            theNeuralEngine.neuralPipeline.coneMosaic.visualize('activation', squeeze(diffResponse), 'verticalActivationColorBarInside', true);
        
            % Also visualize the full absorptions density
            figNo = 999;
            theNeuralEngine.neuralPipeline.coneMosaic.visualizeFullAbsorptionsDensity(figNo);
        end
        
        if (isfield(theNeuralEngine.neuralPipeline, 'mRGCmosaic'))
            theNeuralEngine.neuralPipeline.mRGCmosaic.visualizeResponses(...
                responseTemporalSupportSeconds, inSampleTestStimResponses(trainNoiseFlag), ...
                'stimulusTemporalSupportSeconds', temporalSupport,...
                'stimulusSceneSequence', testScene);
        end
        
    end
    
    % Train the classifier. This shows the usage to extact information
    % from the container retrned as the first return value from the neural
    % response engine - we index the responses by the string contained in
    % the variable trainFlag (which was itself passed to the neural
    % repsonse engine above.)
    %
    % Once extracted from the container, the responses are a 3 dimensional
    % matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
    %   instancesNum   - number of response instances
    %   mNeuralDim     - dimension of neural response at one timepoint
    %   tTimeBins      - number of time points in stimulus sequence.
    theClassifierEngine.compute('train', ...
        inSampleNullStimResponses(trainNoiseFlag), ...
        inSampleTestStimResponses(trainNoiseFlag));
    
    % Save computed response instances
    if (saveResponses)
        responses.inSampleNullStimResponses = inSampleNullStimResponses;
        responses.inSampleTestStimResponses = inSampleTestStimResponses;
    end
    
end

% Predict using trained classifier.
%
% Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% indicates what type of noise is used to generate the stimuli used for
% prediction.  Typically 'random'.
[outOfSampleNullStimResponses, ~] = theNeuralEngine.compute(...
    nullScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testNoiseFlag});

% Generate stimuli for prediction, TEST stimulus
[outOfSampleTestStimResponses, ~] = theNeuralEngine.compute(...
    testScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testNoiseFlag});

% Do the prediction
dataOut = theClassifierEngine.compute('predict', ...
    outOfSampleNullStimResponses(testNoiseFlag), ...
    outOfSampleTestStimResponses(testNoiseFlag));

% Save computed response instances
if (saveResponses)
    responses.outOfSampleNullStimResponses = outOfSampleNullStimResponses;
    responses.outOfSampleTestStimResponses = outOfSampleTestStimResponses;
end
    
% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end