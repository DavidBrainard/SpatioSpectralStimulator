function [modulatingPrimary,upperPrimary,lowerPrimary] = ReceptorIsolate(T_receptors,whichReceptorsToIsolate, ...
    whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary, ...
    primaryHeadroom,maxPowerDiff,targetBasis,projectIndices,targetLambda,desiredContrasts,ambientSpd)
% [modulatingPrimary,upperPrimary,lowerPrimary] = ReceptorIsolate(T_receptors,whichReceptorsToIsolate, ...
%   whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary, ...
%   primaryHeadroom,maxPowerDiff,targetBasis,projectIndices,targetLambda,[desiredContrasts],[ambientSpd])
%
% Find the best isolating modulation around a given background.  This is a very general routine,
% with inputs as follows.
%
% T_receptors -             Spectral sensitivities of all receptors being considered, in standard PTB format.
% whichReceptorsToIsolate - Index vector specifying which receptors we want to modulate.
% whichReceptorsToIgnore -  Index vector specifying receptors where we don't care what they do. Can be the empty matrix.
%                           Why, you ask, might you want to do this?  Maybe if T_receptors contains the rods but you're
%                           working at high light levels.
% B_primary -               These calculations are device dependent.  B_primary is a set of basis vectors for the lights
%                           that the device can produce, scaled so that the gamut is for the range [0-1] on each primary.
% backgroundPrimary -       Background around which modulation will occur, in primary space.
% initialPrimary -          Initial guess for primary search.  Usually passed as backgroundPrimary.
% primaryHeadroom -         If you don't want to modulate all the way to the edge of the device gamut, pass a number in the
%                           range [0-1] here.  This constrains the primary settings to be within [0+primaryHeadRoom,1-primaryHeadRoom].
%                           This can be useful for anticipating the fact that the device may get dimmer over time, so that the
%                           stimulus you compute at time 0 remains in gamut after a subsequent calibration.
% maxPowerDiff -            This enforces a smoothness constraint on the spectrum of the computed modulation.  You wouldn't
%                           use this for a device like a monitor, but for our OneLight device this prevents solutions that
%                           wiggle rapdily as a function of wavelength.  Our intuition is that such wiggly spectra are not
%                           as robust in their isolationg properties as smooth spectra.  Pass Inf to ignore.
% targetBasis -             Keep produced spectrum close to the subspaced defined by this basis set.
% projectIndices -          Only wavelengths corresponding to these indices count in the subspace evaluation.  Pass empty matrix to use them all.
% targetLambda -            In trying to match subspace, weight the targetSpectrum difference by this scalar. Set to zero to ignore subspace.
% desiredContrasts -        Vector of target contrasts for receptors that will be isolated.  This is useful, for example,
%                           if you want to do something like produce a modulation with equal L and M cone contrasts with
%                           opposite signs while at the same time silencing the S cones.  This vector should have the same
%                           length as whichReceptorsToIsolate.  It can be the empty vector, in which case the routine maximizes
%                           the sum of the contrasts of the receptors in whichReceptorsToIsolate.
% ambientSpd -              Spectral power distribution of the ambient light.  Optional.  Defaults to zero.
%
% Contrast held at zero for any receptor classes not in the isolate/ignore lists.
%
% Known Bugs:
%   A) It looks like the code that enforces gamut limitations in a manner that
%   handles backgrounds that do not correspond to device primary settings
%   of 0.5 only works just right if no primaries are being pinned.  There
%   is an error check at the end of the function which throws an error if
%   any of the primary values returned are outside the range [0-1], so our
%   motivation for thinking about this will only cross threshold if this
%   error ever gets thrown.
%
% 05/15/20   dhb      Wrote it from ReceptorIsolate

%% Check and if necessary set default for projectIndices
if (isempty(projectIndices))
    projectIndices = 1:size(targetBasis,1);
end

%% Check whether the desired contrasts were passed, and if so check
% consistency of its dimensions.
if (nargin < 11)
    desiredContrasts = [];
end
if ~isempty(desiredContrasts)
    if length(whichReceptorsToIsolate) ~= length(desiredContrasts)
        error('Size of whichReceptorsToIsolate and of desired contrasts vector do not line up')
    end
end

%% Default for ambientSpd
if (nargin < 12 || isempty(ambientSpd))
    ambientSpd = zeros(size(B_primary,1),1);
end

%% Initial guess for modulation
x = initialPrimary-backgroundPrimary;

%% Figure out which receptors get zero modulation and set up constraint for this.
whichReceptorsToZero = setdiff(1:size(T_receptors,1),[whichReceptorsToIsolate whichReceptorsToIgnore]);
backgroundReceptors = T_receptors*B_primary*backgroundPrimary;
backgroundReceptorsZero = backgroundReceptors(whichReceptorsToZero);
Aeq = T_receptors(whichReceptorsToZero,:)*B_primary;
beq = backgroundReceptorsZero;

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

C = [eye(size(initialPrimary,1)) ; ...
    -eye(size(initialPrimary,1)) ; ...
    -eye(size(initialPrimary,1)) ; ...
    eye(size(initialPrimary,1))];
q = [ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
     zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ; ...
     ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
     zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ];
 
% C = [eye(size(initialPrimary,1)) ; ...
%     -eye(size(initialPrimary,1)) ; ...
%     ];
% q = [ones(length(initialPrimary),1) - primaryHeadroom - backgroundPrimary ; ...
%      zeros(length(initialPrimary),1) - primaryHeadroom + backgroundPrimary ; ...
%      ];


%% Fix numerical issues with vlb > vub that can sometimes come up.
vlbTolerance = 1e-6;
for ii = 1:length(vub)
    if (vlb(ii) > vub(ii) - vlbTolerance)
        vlb(ii) = vub(ii) - vlbTolerance;
    end
end

%% Optimize.
% Progressive smoothing seems to work better than providing final value all
% at once.
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','on','Algorithm','sqp', 'MaxFunEvals', 100000, 'TolFun', 1e-10, 'TolCon', 1e-10, 'TolX', 1e-5);

modulatingPrimary = fmincon(@(x) IsolateFunction(x,B_primary,backgroundPrimary,ambientSpd,T_receptors,whichReceptorsToIsolate,desiredContrasts,targetBasis,projectIndices,targetLambda), ...
    x,C,q,Aeq,beq,[],[],[],options);

%% Check gamut
upperPrimary = backgroundPrimary + modulatingPrimary;
lowerPrimary = backgroundPrimary - modulatingPrimary;
primaryTolerance = 1e-6;
if (any(upperPrimary > 1-primaryHeadroom+primaryTolerance))
    error('upperPrimary too large');
end
if (any(upperPrimary < 0+primaryHeadroom-primaryTolerance))
    error('upperPrimary too small');
end


end

% f = IsolateFunction(x,B_primary,backgroundPrimary,T_receptors,whichReceptorsToIsolate,C,lambda)
%
% Optimization subfunction.  This mixes maximizing response of isolated
% receptors with smoothness.
function f = IsolateFunction(x,B_primary,backgroundPrimary,ambientSpd,T_receptors,whichReceptorsToIsolate,desiredContrasts,...
    targetBasis,projectIndices,targetLambda)

% Compute background including ambient
backgroundSpd = B_primary*backgroundPrimary + ambientSpd;
modulationSpd = B_primary*x + backgroundSpd;
isolateContrast = ExcitationsToContrast(T_receptors(whichReceptorsToIsolate,:)*modulationSpd, ...
    T_receptors(whichReceptorsToIsolate,:)*backgroundSpd);

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
    f1 = sqrt( sum((isolateContrast-desiredContrasts').^2) ) / norm(desiredContrasts);
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