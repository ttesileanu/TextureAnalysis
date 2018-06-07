function [gain, predictions, details] = getTernaryPredictions(...
    ni_ev, pp_data, covIn, covOut, lagrange, varargin)
% getTernaryPredictions Calculate predictions for the ternary texture case
% starting from natural image patches.
%   [gain, predictions] = ...
%       getTernaryPredictions(ni_ev, pp_data, covIn, covOut, lagrange)
%   returns the `gain` matrix obtained from an efficient coding argument on
%   the natural image statistics `ni_ev`, assuming covariance matrices for
%   input and output noise `covIn` and `covOut` respectively, and a
%   Lagrange multiplier `lagrange` (see solveLinearEfficientCoding for
%   details). It also uses this gain matrix to generate threshold
%   predictions in the directions indicated by the `groups` and
%   `directions` fields of the structure `pp_data`, and finally scales
%   these to optimally match the `thresholds` from `pp_data` (see fitscale
%   for details).
%
%   [..., details] = getTernaryPredictions(...) returns some details of the
%   calculation, including the `a` and `mse` outputs of fitscale.
%
%   getTernaryPredictions('bootstrap', n) uses `n` bootstrap samples to
%   calculate threshold predictions together with error bars. The
%   arithmetic mean of the bootstrapped gain matrices is returned for
%   `gain`, and the median of the bootstrapped predictions is returned for
%   `predictions`, after using `fitscale` on it. The full list of
%   bootstrapped values is returned in the `details` structure.
%
%   Options:
%    'bootstrap'
%       The number of bootstrap samples to use to estimate predicted
%       thresholds and their error bars. Set to 1 to do no bootstrapping
%       and use the entire dataset.
%    'fitscale_opts'
%       Options to be passed to fitscale.
%    'effcode_opts'
%       Options to be passed to solveLinearEfficientCoding.
%
%   See also: solveLinearEfficientCoding, gainsToThresholds, fitscale.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('bootstrap', 20, @(n) isnumeric(n) && isscalar(n) && n >= 1);
parser.addParameter('fitscale_opts', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('effcode_opts', {}, @(c) iscell(c) && isvector(c));

% 'exclude', strcmp(ternary_avg.groups, 'A_1') & cellfun(@length, ternary_avg.groups) > 6, ...

if nargin == 1 && strcmp(pred, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% start the bootstrapping, store the samples
n_patches = size(ni_ev, 1);

gains = cell(params.bootstrap, 1);
unscaled_predictions = cell(params.bootstrap, 1);
ni_covs = cell(params.bootstrap, 1);

% directions extended to the full dimensionality of texture space
directionsExt = cellfun(@(v) v - 1/3, ternaryextdir(pp_data.groups, pp_data.directions), ...
    'uniform', false);

progress = TextProgress('bootstrapping');
for i = 1:params.bootstrap
    % if there's only one sample, use the whole dataset
    if params.bootstrap == 1
        ni_sample = ni_ev;
    else
        ni_sample = ni_ev(randi(n_patches, n_patches, 1), :);
    end
    
    ni_covs{i} = cov(ni_sample);
    
    gains{i} = solveLinearEfficientCoding(ni_covs{i}, covIn, covOut, lagrange, params.effcode_opts{:});
    unscaled_predictions{i} = gainsToThresholds(gains{i}, directionsExt);
    
    progress.update(100*i/params.bootstrap);
end
progress.done;

% find the average gain and average unscaled prediction
gain = mean(cat(3, gains{:}), 3);
unscaled_predictions_3d = cat(3, unscaled_predictions{:});
unscaled_prediction = median(unscaled_predictions_3d, 3);

% now use fitscale to match predictions to measurements
[acoeff, predictions, pred_mse] = fitscale(...
    unscaled_prediction, pp_data.thresholds, 'stds', pp_data.threshold_intervals, ...
    params.fitscale_opts{:});

% calculate error bars if we have enough samples
if params.bootstrap > 2
    unscaled_lo = quantile(unscaled_predictions_3d, 0.159, 3);
    unscaled_hi = quantile(unscaled_predictions_3d, 0.841, 3);
else
    unscaled_lo = nan(size(unscaled_prediction));
    unscaled_hi = nan(size(unscaled_prediction));
end

predictions_lo = scale(unscaled_lo, acoeff);
predictions_hi = scale(unscaled_hi, acoeff);

% prepare the details structure
details.ni_covs = ni_covs;
details.gains = gains;
details.unscaled_predictions = unscaled_predictions;
details.acoeff = acoeff;
details.pred_mse = pred_mse;
details.prediction_intervals = [predictions_lo(:) predictions_hi(:)];

end

function v = scale(v, a)
% SCALE Scale according to fitscale return factor.

if length(a) == 1
    v = a*v;
else
    mask = (v > 0);
    v(mask) = a(1)*v(mask).^a(2);
    v(~mask) = nan;
end

end