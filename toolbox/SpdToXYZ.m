function [XYZ, xy] = SpdToXYZ(spd,options)
% Calculate CIE XYZ tri-stimulus values from SPDs.
%
% Syntax:
%    [XYZ, xy] = SpdToXYZ(spd)
%
% Description:
%     Calculate CIE XYZ values from the SPDs. It also calculates the CIE
%     (x,y) chromaticity too and plot it with standard color gamuts, sRGB
%     and Adobe RGB.
%
% Inputs:
%     spd -                            Target SPDs to calculate XYZ values.
%
% Outputs:
%     XYZ -                            Absolute XYZ values where Y matches
%                                      to luminance (cd/m2).
%     xy -                             CIE (x, y) chromaticity.
%
% Optional key/value pairs:
%    S -                               Spectrum wavelengths range that the
%                                      spds was measured. Default to
%                                      [380 2 201].
%    verbose -                         Boolean. Default true. Controls
%                                      printout.
%
% See also
%   N/A

% History:
%   02/16/22  smo                      Started On it.

%% Set parameters.
arguments
    spd
    options.S (1,3) = [380 2 201]
    options.verbose (1,1) = true
end

%% Calculate XYZ here.
%
% Load the CMF data and match the spectrum wavelength range.
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos, 683 * T_xyzJuddVos, options.S);

% Calculate it.
XYZ = T_xyz * spd;

%% Calculate CIE xy chromaticity.
xyY = XYZToxyY(XYZ);
xy = xyY(1:2,:);

%% Plot it if you want.
%
% Some standard color gamuts and spectral locus to compare.
% sRGB gamut.
sRGBGamut = [0.6400 0.3000 0.1500; 0.3300 0.6000 0.0600];
sRGBGamut(:,end+1) = sRGBGamut(:,1);

% Adobe RGB gamut.
adobeRGBGamut = [0.6400, 0.2100, 0.1500; 0.3300 0.7100 0.0600];
adobeRGBGamut(:,end+1) = adobeRGBGamut(:,1);

% Set Spectral locus.
spectralLocus = XYZToxyY(T_xyz);
spectralLocus(:,end+1) = spectralLocus(:,1);

% Plot here.
if (options.verbose)
    figure; clf; hold on;
    plot(xyY(1,:), xyY(2,:), 'r.', 'MarkerSize', 10);
    plot(sRGBGamut(1,:), sRGBGamut(2,:), 'g-', 'LineWidth', 1);
    plot(adobeRGBGamut(1,:), adobeRGBGamut(2,:), 'b-', 'LineWidth', 1);
    plot(spectralLocus(1,:), spectralLocus(2,:), 'k-');
    xlabel('CIE x');
    ylabel('CIE y');
    legend('Test', 'sRGB', 'AdobeRGB', 'Spectral locus');
end

end