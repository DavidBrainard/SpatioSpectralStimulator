%% Generate a custom experimental display
    % Optional settable key-value pairs are:
    % 'dotsPerInch'                                             - scalar
    % 'viewingDistanceMeters'                                   - scalar
    % 'gammaTable'                                              - [mValues x 3] matrix of LUTs
    % 'wavelengthSupportNanoMeters'                             - [nWaves x 1] matrix of wavelengths
    % 'ambientSPDWattsPerSteradianM2NanoMeter'                  - [nWaves x 1] matrix of the ambient SPD
    % 'spectralPowerDistributionWattsPerSteradianM2NanoMeter'   - [nWaves x 3] matrix of the RGB guns SPDs
    
    presentationDisplay = generateCustomDisplay(...
        'dotsPerInch', 220, ...
        'viewingDistanceMeters', 1.00, ...
        'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
        'plotCharacteristics', true);