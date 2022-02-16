function [XYZ] = SpdToXYZ(spd,options)
% Calculate CIE XYZ tri-stimulus values from SPDs.
%
% Syntax:
%    [XYZ] = SpdToXYZ(spd)
%
% Description:
%     Calculate CIE XYZ values from the SPDs.
%
% Inputs:
%     spd -                            Target SPDs to calculate XYZ values.
%
% Outputs:
%     XYZ -                            Absolute XYZ values where Y matches
%                                      to luminance (cd/m2). 
%
% Optional key/value pairs:
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
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos, 683 * T_xyzJuddVos, options.S);

% Calculate XYZ values.
XYZ = T_xyz * spd;

% CIE xy chromaticity. It can be printed out too if you want. Not now.
xyY = XYZToxyY(XYZ);

end