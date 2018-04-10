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
%
%   Options:
%    'exclude'
%       Mask for positions in the `pred` and `meas` vectors to ignore when
%       finding the optimal factor `a`.
%    'stds'
%       Error bars for each measurement.
%    'log'
%       If true, the fit instead attempts to optimize the difference
%       between `log(predscaled)` and `log(meas)`, where
%           predscaled = a*pred .
%       In this case, the errors are assumed to be distributed according to
%       a log-normal, so that the objective function that is to be
%       minimized (and that occurs in the return argument `mse`) is
%           mse = mean( ((log(predscaled) - log(meas)) ./ stds) .^2 )
%    'logslope'
%       If true, the fit uses both a slope and an intercept in log space,
%           predscaled = a(1)*pred.^a(2) .

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

if nargin == 1 && strcmp(pred, 'defaults')
    show_defaults = true;
    meas = [];
else
    show_defaults = false;
end

parser.addParameter('stds', [], @(v) isnumeric(v) && isvector(v) && length(v) == length(meas));
parser.addParameter('exclude', false(size(meas)), @(b) islogical(b) && isvector(b) && length(b) == length(meas));
parser.addParameter('log', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('logslope', false, @(b) islogical(b) && isscalar(b));

if show_defaults
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

if params.log
    % need the values to also be positive
    meas_mask = meas_mask & (meas > 0);
    pred_mask = pred_mask & (pred > 0);
end

mask = meas_mask & pred_mask;

if ~params.log
    % dividing by error bars maps to problem where error bars are constant
    prednorm = pred(mask) ./ params.stds(mask);
    measnorm = meas(mask) ./ params.stds(mask);

    % find scaling factor
    a = dot(prednorm, measnorm) / dot(prednorm, prednorm);
    predscaled = a*pred;
else
    x = log(pred(mask));
    y = log(meas(mask));
    logweights = 1./params.stds(mask).^2;
    if ~params.logslope
        loga = wmean(y - x, logweights);
        xscaled = x + loga;
        a = exp(loga);
        predscaled = a*pred;
    else
        xwm = wmean(x, logweights);
        ywm = wmean(y, logweights);
        xywm = wmean(x.*y, logweights);
        x2wm = wmean(x.^2, logweights);
        
        a(2) = (xywm - xwm*ywm) / (x2wm - xwm^2);
        loga = wmean(y - a(2)*x, logweights);
        
        xscaled = a(2)*x + loga;
        
        a(1) = exp(loga);
        
        pos_mask = (pred > 0);
        predscaled = nan(size(pred));
        predscaled(pos_mask) = a(1)*pred(pos_mask).^a(2);
    end
end

if nargout > 2
    if sum(mask) > 0
        % if we don't have predictions for any of the entries where we have
        % good measurements, we're doing a bad fit
        if sum(meas_mask & ~pred_mask) > 0
            mse = inf;
        else
            if ~params.log
                mse = mean( ((predscaled(mask) - meas(mask)) ./ params.stds(mask)) .^2 );
            else
                mse = wmean((xscaled - y) .^ 2, logweights);
            end
        end
    else
        mse = nan;
    end
end

end

function m = wmean(v, w)
% WMEAN Weighted mean.

m = mean(v .* w) / mean(w);

end