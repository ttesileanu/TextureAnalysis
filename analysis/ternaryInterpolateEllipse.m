function [mags_int, uvecs_int] = ternaryInterpolateEllipse(mags, uvecs, n)
% ternaryInterpolateEllipse Interpolate threshold locations in a ternary
% texture group assuming an elliptic threshold surface.
%   [mags_int, uvecs_int] = ternaryInterpolateEllipse(mags, uvecs, n)
%   generates `n` interpolated thresholds along the ellipse that best fits
%   the given thresholds `mags` in the directions `uvecs`. The fit is
%   performed after projecting to the 2-component space (see TERNARY3TO2).
%
%   If a fit isn't possible because there are not enough finite data
%   points, or if the fit isn't an ellipse, the function returns empty
%   matrices.
%
%   The input `uvecs` can be a cell array or 3-component vectors or a
%   3-column matrix. The returned `uvecs_int` is always a matrix.

% convert input vectors to matrix
if iscell(uvecs)
    uvecs = cell2mat(cellfun(@(v) v(:)', uvecs(:), 'uniform', false));
end

% make sure we focus on valid inputs
mask = isfinite(mags) & all(isfinite(uvecs), 2);

if sum(mask) < 3
    % not enough points for an ellipse fit
    mags_int = [];
    uvecs_int = [];
    return;
end

mags = mags(mask);
uvecs = uvecs(mask, :);

% project to 2d
projected = ternary3to2(ternaryrec(mags, uvecs));

% try to fit ellipse
M = fit_ellipse(projected);
tol = 1e-6;
if min(eig(M)) < -tol
    % not positive definite, not an ellipse
    mags_int = [];
    uvecs_int = [];
    return;
end

% generate directions
phi = linspace(0, 2*pi, n+1);
phi = phi(1:end-1); % 2*pi is the same as 0
v = [cos(phi) ; sin(phi)];

% generate new projected threshold locations
m = sqrt(1 ./ sum(v.* (M*v), 1));
projected_int = bsxfun(@times, m, v);

% project back to quasi-probabilities
[mags_int, uvecs_int] = ternarydec(ternary2to3(projected_int'));

end