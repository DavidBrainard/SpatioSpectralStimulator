function [peakWavelengths] = FindPeakSpds(spds,options)
% Find a peak wavelength within the given spds.
%
% Syntax:
%    [peakWavelengths] = FindPeakSpds(spds)
%
% Description:
%    This finds a peak wavelength for the given spectral power
%    distributions. It is useful to use finding a peak wavelength for LEDs.
%    Also, we can utilize this when we use power meter to set its measuring
%    wavelength.
%
% Inputs:
%    spds                       - Spectral power distrubutions that you
%                                 want to find its peak wavelength. You can
%                                 pass this as array to find multiple peaks
%                                 for the spds more than two.
%
% Outputs:
%    peakWavelengths            - Peak wavelength for the given spds.
%
% Optional key/value pairs:
%    S                          - Default to [380 2 201]. The wavelength
%                                 range of the test spectrum.
%    verbose                    - Default to rue. Boolean. Controls the
%                                 printout.
%
% See also:
%    N/A

% History:
%    03/21/22  smo              - Started on it

%% Set parameters.
arguments
    spds
    options.verbose (1,1) = true
    options.S (1,3) = [380 2 201]
end

%% Find a peak wavelength here.
%
% Set the number of spds to test.
nSpds = size(spds,2);
spdStart = options.S(1);
spdInterval = options.S(2);

% Make a loop to find a peak wavelength per each test spectrum.
for ss = 1:nSpds
    % Find a peak wavelength here and save the results in an array. The
    % results are available in peakWavelegnths.
    spdSingle = spds(:,ss);
    peakWavelength = max(spdSingle);
    
    % If we found more than one wavelengths that have the peak, we take the first one.
    idxPeakWavelength(ss) = min(find(spdSingle(:) == peakWavelength));
    peakWavelengths(ss) = spdStart + spdInterval * (idxPeakWavelength(ss) - 1);
    
    % Print out the results if you want.
    if (options.verbose)
        fprintf('Spectrum (%2.0d) peak wavelength: %d nm \n', ss, peakWavelengths(ss));
    end
end

% Plot the results if you want.
if (options.verbose)
    wls = SToWls(options.S);
    
    figure; clf; hold on;
    % Plot the all test spectra.
    plot(wls,spds);
    
    % Mark peak wavelength and also make legend in cell array to add in the
    % plot.
    for ss = 1:nSpds
        plot(peakWavelengths(ss),spds(idxPeakWavelength(ss),ss),'o',...
            'MarkerSize',8,'MarkerEdgeColor',zeros(3,1),'MarkerFaceColor',[1 0 0]);
        legends{ss} = append('Spd',num2str(ss));
    end
    legends{end+1} = 'Peak';
    xlabel('Wavelength (nm)','FontSize',15);
    ylabel('Spectral irrdiance','FontSize',15);
    xlim([380 780]);
    legend(legends','FontSize',10);
end

end