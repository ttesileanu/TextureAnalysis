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
%   Directions for which either the predicted value or its error bar (see
%   'stds' option below) are infinite or NaN are ignored. The default
%   metric used for comparing the thresholds is mean of squared differences.
%
%   [..., meta_details, details] = optimizePredictedThresholds(...) also
%   returns a structure of convergence details for the `radius2`
%   optimization, and a cell array of convergence results as returned from
%   `mapIntersectEllipsoid`.
%
%   Options:
%    'exclude'
%       Boolean mask of directions that should be ignored for the fit.
%    'stds'
%       Uncertainty values for the measurements. Measurements with infinite
%       or NaN error bars are ignored for the fit. When uncertainty values
%       are present, the differences between predictions and measurements
%       are divided by these uncertainties before being used in the loss
%       function. This can also be given by an N x 2 matrix where the first
%       column is the lower bound of the error interval, and the second
%       column is the upper bound.
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
parser.addParameter('exclude', [], @(v) islogical(v) && isvector(v));
parser.addParameter('stds', [], @(v) isnumeric(v) && ((isvector(v) && length(v) == length(interpolation)) ...
    || (ismatrix(v) && size(v, 2) == 2 && size(v, 1) == length(interpolation))));
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

% handle default exclude mask
if isempty(params.exclude)
    params.exclude = false(size(interpolation));
end

% default to constant error bars
if isempty(params.stds)
    params.stds = ones(size(interpolation));
elseif ~isvector(params.stds)
    params.stds = diff(params.stds, [], 2);
end

% mask of values to use
use_mask = (~params.exclude(:) & isfinite(expected(:)));

% calculate the MSE between the two vectors, normalizing by error bars, and
% ignoring entries that are infinite or NaN
mse_mask = @(pred, meas, s, mask) mean( ((pred(mask) - meas(mask)) ./ s(mask)).^2 );
mse_fct = @(pred) mse_mask(pred, expected, params.stds, use_mask & isfinite(pred(:)));

% run the optimization
[optimal_radius2, optimal_fval, optimal_exitflag, optimal_output] = fminbnd(@(radius2) ...
    mse_fct(mapIntersectEllipsoid(interpolation, trafo, radius2)), ...
    radius2_range(1), radius2_range(2), params.fminbnd_opts{:});

% recalculate the predicted thresholds at the optimal noise radius
[predictions, mapped, details] = mapIntersectEllipsoid(interpolation, trafo, ...
    optimal_radius2);

meta_details.radius2 = optimal_radius2;
meta_details.min_mse = optimal_fval;
meta_details.exitflag = optimal_exitflag;
meta_details.output = optimal_output;

end