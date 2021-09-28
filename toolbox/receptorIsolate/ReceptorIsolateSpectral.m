function [modulatingPrimary,upperPrimary,lowerPrimary] = ReceptorIsolateSpectral(T_receptors,desiredContrasts, ...
    B_primary,backgroundPrimary,initialPrimary, ...
    primaryHeadroom,targetBasis,projectIndices,targetLambda,ambientSpd,varargin)
% Find modulating primaries to produce desired contrast/excitations re background.
% 
% Syntax:
%    [modulatingPrimary,upperPrimary,lowerPrimary] = ReceptorIsolateSpectral(T_receptors,desiredContrats, ...
%       B_primary,backgroundPrimary,initialPrimary, ...
%       primaryHeadroom,targetBasis,projectIndices,targetLambda,ambientSpd)
%
% Description:
%     Find the best isolating modulation around a given background.  This
%     is a very general routine. It enforces constraints on primary values
%     (between 0 and 1, with possible headroom specified) and can enforce
%     a constraint that maximizes projection onto a subspace speicifed by
%     a set of basis vectors.
%
% Inputs:
%     T_receptors -             Spectral sensitivities of all receptors being considered, in standard PTB format.
%     desiredContrasts -        Vector of target contrasts for receptors that will be isolated.
%     B_primary -               These calculations are device dependent.  B_primary is a set of basis vectors for the lights
%                               that the device can produce, scaled so that the gamut is for the range [0-1] on each primary.
%     backgroundPrimary -       Background around which modulation will occur, in primary space.
%     initialPrimary -          Initial guess for primary search.  Usually passed as backgroundPrimary.
%     primaryHeadroom -         If you don't want to modulate all the way to the edge of the device gamut, pass a number in the
%                               range [0-1] here.  This constrains the primary settings to be within [0+primaryHeadRoom,1-primaryHeadRoom].
%                               This can be useful for anticipating the fact that the device may get dimmer over time, so that the
%                               stimulus you compute at time 0 remains in gamut after a subsequent calibration.
%     targetBasis -             Keep produced spectrum close to the subspaced defined by this basis set.
%     projectIndices -          Only wavelengths corresponding to these indices count in the subspace evaluation.  Pass empty matrix to use them all.
%     targetLambda -            In trying to match subspace, weight the targetSpectrum difference by this scalar. Set to zero to ignore subspace.
%     ambientSpd -              Spectral power distribution of the ambient light.
%
% Outputs:
%    modulatingPrimary -        Primary values to add to background to produce desired contrast/excitations
%    upperPrimary -             The primary values when modulatingPrimary is added to the background.
%    lowerPrimary -             The primary values when modulatingPrimary is subtracted from the background
%
% Optional key/value pairs:
%    POSITIVE_ONLY         - Boolean (default false). Only consider gamut constraints for positive arm of modulation.
%    EXCITATIONS           - Boolean (default false). Aim for target excitations, rather than contrasts. By default it does contrast.
%
% 05/15/20   dhb            Wrote it from ReceptorIsolate
% 09/28/21   dhb, smo       Improve comments.

%% Parse key/value pairs
p = inputParser;
p.addParameter('POSITIVE_ONLY', false, @islogical);
p.addParameter('EXCITATIONS', false, @islogical);
p.parse(varargin{:});

%% Check and if necessary set default for projectIndices
if (isempty(projectIndices))
    projectIndices = 1:size(targetBasis,1);
end

%% Check consitency of passed desired contrasts with photoreceptors
    if (size(T_receptors,1) ~= length(desiredContrasts))
        error('Size of whichReceptorsToIsolate and of desired contrasts vector do not line up')
    end

%% Initial guess for modulation
x = initialPrimary-backgroundPrimary;

% Since our modulations are symmetric, we need to make sure that we're not
% out of gamut if our background is not constant across wl band. For a
% half-on background, both the positive and the negative poles of a
% modulation will be in gamut, but that's not necessary the case if the
% background is not 0.5 for all wl bands.
%
% The following piece of code may also only work just right if we're
% not pinning primaries.
primaryHeadRoomTolerance = 1e-7;
if (any(backgroundPrimary < primaryHeadroom - primaryHeadRoomTolerance))
    error('Cannot work if background primary is less than specified headroom');
end
if (any(backgroundPrimary > 1-primaryHeadroom+primaryHeadRoomTolerance))
    error('Cannot work if background primary is greater than 1 minus specified headroom');
end
for b = 1:size(backgroundPrimary, 1)
    if backgroundPrimary(b) > 0.5
        vub(b) = 1-primaryHeadroom;
        vlb(b) = backgroundPrimary(b)-(1-primaryHeadroom-backgroundPrimary(b));
    elseif backgroundPrimary(b) < 0.5
        vub(b) = backgroundPrimary(b)+(backgroundPrimary(b)-primaryHeadroom);
        vlb(b) = primaryHeadroom;
    elseif backgroundPrimary(b) == 0.5
        vub(b) = 1-primaryHeadroom;
        vlb(b) = primaryHeadroom;
    end
end

POSITIVE_ONLY = p.Results.POSITIVE_ONLY;
if (POSITIVE_ONLY)
    C = [eye(size(initialPrimary,1)) ; ...
        -eye(size(initialPrimary,1)) ; ...
        ];
    q = [ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
        zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ; ...
        ]; 
else
    C = [eye(size(initialPrimary,1)) ; ...
        -eye(size(initialPrimary,1)) ; ...
        -eye(size(initialPrimary,1)) ; ...
        eye(size(initialPrimary,1))];
    q = [ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
        zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ; ...
        ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
        zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ];
end

%% Optimize.
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off'); %,'LargeScale','on','Algorithm','sqp','MaxFunEvals', 100000, 'TolFun', 1e-10, 'TolCon', 1e-10, 'TolX', 1e-5);

modulatingPrimary = fmincon(@(x) IsolateFunction(x,B_primary,backgroundPrimary,ambientSpd, ...
    T_receptors,desiredContrasts,targetBasis,projectIndices,targetLambda,p.Results.EXCITATIONS), ...
    x,C,q,[],[],[],[],[],options);

%% Check gamut
upperPrimary = backgroundPrimary + modulatingPrimary;
lowerPrimary = backgroundPrimary - modulatingPrimary;
primaryTolerance = 2e-2;
if (any(upperPrimary > 1-primaryHeadroom+primaryTolerance))
    error('upperPrimary too large');
end
if (any(upperPrimary < 0+primaryHeadroom-primaryTolerance))
    error('upperPrimary too small');
end
if (~POSITIVE_ONLY)
    if (any(lowerPrimary > 1-primaryHeadroom+primaryTolerance))
        error('lowerPrimary too large');
    end
    if (any(lowerPrimary < 0+primaryHeadroom-primaryTolerance))
        error('lowerPrimary too small');
    end
end

end

% f = IsolateFunction(x,B_primary,backgroundPrimary,ambientSpd,T_receptors,desiredContrasts,...
%         targetBasis,projectIndices,targetLambda,EXCITATIONS)
%
% Optimization subfunction.  This mixes maximizing acheiving
% desiredContrasts with staying within the target basis.
function f = IsolateFunction(x,B_primary,backgroundPrimary,ambientSpd,T_receptors,desiredContrasts,...
    targetBasis,projectIndices,targetLambda,EXCITATIONS)

% Compute background including ambient
backgroundSpd = B_primary*backgroundPrimary + ambientSpd;
modulationSpd = B_primary*x + backgroundSpd;
if (EXCITATIONS)
    isolateContrast = T_receptors*modulationSpd;
else
    isolateContrast = ExcitationsToContrast(T_receptors*modulationSpd,T_receptors*backgroundSpd);
end

% Get contrast term
if isempty(desiredContrasts)
    % Want the sum of the isolated receptor contrasts to be big. fmincon
    % minimizes, hence the negative sign.  Acheive this by minimizing
    % the difference between the isolatedContrasts and unity.  For
    % reasons not fully understood, this works better numerically than
    % simply minimizing the negative sum of squared contrasts.
    f1 = sqrt( sum((isolateContrast-1).^2) ) / norm(desiredContrasts);
else
    % Minimize difference between desired and what we get
    f1 = sqrt( sum((isolateContrast-desiredContrasts).^2) ) / norm(desiredContrasts);
end

% Get basis term
if (targetLambda > 0)
    projectionSpd = targetBasis*(targetBasis(projectIndices,:)\modulationSpd(projectIndices));
    f2 = sqrt(sum( (modulationSpd(projectIndices)-projectionSpd(projectIndices)).^2 ) ) / norm(modulationSpd(projectIndices));
    f = f1 + targetLambda*f2;
else
    f = f1;
end

end