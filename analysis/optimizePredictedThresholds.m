function [predictions, mapped, meta_details, details] = ...
    optimizePredictedThresholds(expected, interpolation, trafo, radius2_range, varargin)
% optimizePredictedThresholds Optimize the noise radius to find the best
% match between predicted and measured thresholds.
%   [predictions, mapped] = optimizePredictedThresholds(...
%       expected, interpolation, trafo, [min_radius2, max_radius2])
%   searches for the `radius2` parameter in `mapIntersectEllipsoid` in the
%   given `[min_radius2, max_radius2]` range until it obtains a best match
%   between the absolute values of the predicted thresholds (the
%   intersection points obtained from running `mapIntersectEllipsoid` on
%   `interpolation`, given the matrix `trafo`) and the `expected` values.
%   The predicted thresholds at the optimum are returned in `predictions`;
%   their mappings through `interpolation` are returned in `mapped`.
%   Directions for which either the `expected` value or the predicted one
%   are infinite or NaN are ignored. The default metric used for comparing
%   the thresholds is sum of squared differences.
%
%   [..., meta_details, details] = optimizePredictedThresholds(...) also
%   returns a structure of convergence details for the `radius2`
%   optimization, and a cell array of convergence results as returned from
%   `mapIntersectEllipsoid`.
%
%   Options:
%    'ignore_mask'
%       Boolean mask that should contain `true` for directions that should
%       be ignored for the fit. The predictions for these directions will
%       still be returend by the function.
%    'fminbnd_opts'
%       Options to pass to FMINBND.
%
%   See also: mapIntersectEllipsoid, FMINBND.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% generate default fsolve options
default_opts = optimset('display', 'iter', 'tolx', 1e-10);

% add possible parameters with their defaults
parser.addParameter('ignore_mask', [], @(v) islogical(v) && isvector(v));
parser.addParameter('fminbnd_opts', {default_opts}, @(c) iscell(c));

% handle showing of defaults.
if nargin == 1 && strcmp(expected, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse optional arguments
parser.parse(varargin{:});
params = parser.Results;

% handle default ignore mask
if isempty(params.ignore_mask)
    params.ignore_mask = false(size(interpolation));
end

% mask of values to use
use_mask = (~params.ignore_mask(:) & isfinite(expected(:)));

% calculate the L2 norm of the difference between two vectors while
% ignoring entries that are infinite or NaN in either of them
normdiff2_mask = @(a, b, mask) norm((a(mask) - b(mask)) / (sum(mask) - 1))^2;
normdiff2_skipinf = @(a, b) normdiff2_mask(a(:), b(:), use_mask(:) & isfinite(a(:)) & isfinite(b(:)));
normdiff2_skipinf_abs = @(a, b) normdiff2_skipinf(abs(a), abs(b));

% % calculate the Spearman correlation between two vectors while ignoring
% % entries that are infinite or NaN in either of them
% spcorr_mask = @(a, b, mask) corr(a(mask), b(mask), 'type', 'spearman');
% spcorr_skipinf = @(a, b) spcorr_mask(a(:), b(:), ignore_mask & isfinite(a) & isfinite(b));

% run the optimization
[optimal_radius2, optimal_fval, optimal_exitflag, optimal_output] = fminbnd(@(radius2) ...
    normdiff2_skipinf_abs(...
        mapIntersectEllipsoid(interpolation, trafo, radius2), expected), ...
    radius2_range(1), radius2_range(2), params.fminbnd_opts{:});

% recalculate the predicted thresholds at the optimal noise radius

[predictions, mapped, details] = mapIntersectEllipsoid(interpolation, trafo, ...
    optimal_radius2);

meta_details.radius2 = optimal_radius2;
meta_details.min_mse = optimal_fval;
meta_details.exitflag = optimal_exitflag;
meta_details.output = optimal_output;

end