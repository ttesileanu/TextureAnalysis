function [a, predscaled, mse] = fitscale(pred, meas, varargin)
% FITSCALE Find the best scaling factor to match a vector of predictions to
% data.
%   a = FITSCALE(pred, meas) finds the factor `a` such that
%   norm(a*pred - meas) is minimized. Entries whether either measurements
%   or their error bars (see below) are infinite or NaN are ignored. If
%   some entries of the predictions are infinite or NaN when measurements
%   and error bars aren't, these are also ignored, but the mean-squared
%   error (see below) that is returned is infinite.
%
%   [a, predscaled] = FITSCALE(pred, meas) also returns the scaled
%   predictions, predscaled = a*pred.
%
%   FITSCALE(pred, meas, 'stds', stds) assumes that the measured values
%   have error bars as given in the `stds` vector. Measurements with
%   infinite or NaN error bars are ignored.
%
%   FITSCALE(..., 'exclude', mask) uses the binary `mask` to exclude some
%   of the entries from `pred` and `meas` from the fit.
%
%   [..., mse] = FITSCALE(...) also returns the mean squared error
%   obtained with the optimal scaling coefficient. Specifically,
%       mse = mean( ((a*pred - meas) ./ stds) .^2 )
%   If error bars aren't given, `stds` is assumed to be identically 1. If
%   any predictions are infinite or NaN at positions where measurements and
%   their error bars are finite, the `mse` is returned as `inf`.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('stds', [], @(v) isnumeric(v) && isvector(v) && length(v) == length(meas));
parser.addParameter('exclude', false(size(meas)), @(b) islogical(b) && isvector(b) && length(b) == length(meas));

if nargin == 1 && strcmp(pred, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% default to constant error bars
if isempty(params.stds)
    params.stds = ones(size(meas));
end

% figure out which entries to care about
meas_mask = isfinite(meas) & isfinite(params.stds) & ~params.exclude;
pred_mask = isfinite(pred) & ~params.exclude;

mask = meas_mask & pred_mask;

% dividing by error bars maps to problem where error bars are constant
prednorm = pred(mask) ./ params.stds(mask);
measnorm = meas(mask) ./ params.stds(mask);

% find scaling factor
a = dot(prednorm, measnorm) / dot(prednorm, prednorm);
predscaled = a*pred;

if nargout > 2
    if sum(mask) > 0
        % if we don't have predictions for any of the entries where we have
        % good measurements, we're doing a bad fit
        if sum(meas_mask & ~pred_mask) > 0
            mse = inf;
        else
            mse = mean( ((predscaled(mask) - meas(mask)) ./ params.stds(mask)) .^2 );
        end
    else
        mse = nan;
    end
end

end