%% MaxwellianViewCalcs
%
% Refs to Burns and Webb below are to Chapter 28
% of OSA Handbook of Optics, Volume 1.

%% DLP size, converted to horiz and vertical
%
% Use larger dimension field stop size for the
% Maxwellian view system.
chipDiagonalInches = 0.95;
chipDiagonalMm = 25.4*chipDiagonalInches;
chipHToVRatio = 16/9;
chipUnitDiagonal = sqrt(1 + chipHToVRatio^2);
chipHorizMm = chipDiagonalMm/chipUnitDiagonal*chipHToVRatio;
chipVertMm = chipDiagonalMm/chipUnitDiagonal;
check = sqrt(chipHorizMm^2 + chipVertMm^2);
if (abs(chipDiagonalMm - check) > 1e-6)
    error('Did not compute chip size correctly');
end
if (abs(chipHorizMm/chipVertMm - chipHToVRatio) > 1e-6)
        error('Did not compute chip size correctly');
end
stopSizeMm = max([chipHorizMm, chipVertMm]);
fprintf('Field stop size %0.1f mm\n',stopSizeMm);

%% Define optical parameters
%
% Numbering of lenses etc. is outward
% from eye, following Burns & Webb
% Figure 2.
lightDiamSourceMm = 3;

lens1FocalLengthMm = 100;
lens2FocalLengthMm = 100;
eyeFocalLengthMm = 16.67;

lens1DiamMm = 25.4;

%% Compute linear size of field stop on retina in mm
%
% Basically Eqs. 8-9 of Burns and Webb.
retinalSizeMm = (eyeFocalLengthMm/lens1FocalLengthMm)*stopSizeMm;
retinalSizeDeg = 2*atand(retinalSizeMm/(2*eyeFocalLengthMm));
fprintf('Retinal size of stop %0.1f mm, %0.1f deg\n',retinalSizeMm,retinalSizeDeg);

%% Compute diameter of light source in pupil
%
% From Eq. 12 in Burns and Webb.  Makes sense, because
% as the focal length of lens 2 goes to infinity, size
% in pupil goes to zero.
lightDiamPupilMm = (lens1FocalLengthMm/lens2FocalLengthMm)*lightDiamSourceMm;
fprintf('Light source diameter %0.1f mm, light source at pupil %0.1f mm\n',lightDiamSourceMm,lightDiamPupilMm);

%% Compute angular size of lens 1, as seen from center of stop
%
% Assume stop is 1 focal length from lens 1
lens1DiamDegFromStop = 2*atand(lens1DiamMm/(2*lens1FocalLengthMm));
fprintf('Lens 1 diameteter %0.1f mm, %0.1f deg seen from stop\n',lens1DiamMm,lens1DiamDegFromStop);
