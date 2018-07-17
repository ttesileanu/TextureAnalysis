function [a, predScaled, mse] = fitScale(pred, meas, varargin)
% fitScale Find the best scaling factor to match a vector of predictions to
% data.
%   a = fitScale(pred, meas) finds the factor `a` such that
%   `norm(a*pred - meas)` is minimized. Entries whether either measurements
%   or their error bars (see below) are infinite or NaN are ignored. If
%   some entries of the predictions are infinite or NaN when measurements
%   and error bars aren't, these are also ignored, but the mean-squared
%   error (see below) that is returned is infinite.
%
%   [a, predScaled] = fitScale(pred, meas) also returns the scaled
%   predictions, `predScaled = a*pred`.
%
%   fitScale(pred, meas, 'stds', stds) assumes that the measured values
%   have error bars as given in the `stds` vector. This can also be a
%   2-column matrix of measurement intervals. Measurements with infinite or
%   NaN error bars are ignored.
%
%   fitScale(..., 'mask', mask) uses only the measurements and predictions
%   that fit the mask to calculate the scaling coefficient `a`. The
%   returned `predScaled` still contains all the entries, however.
%
%   [..., mse] = fitScale(...) also returns the mean squared error obtained
%   using the optimal scaling coefficient. Specifically,
%       mse = mean( ((a*pred - meas) ./ stds) .^2 )
%   If error bars aren't given, `stds` is assumed to be identically 1. If
%   any predictions are infinite or NaN at positions where measurements and
%   their error bars are finite, the `mse` is returned as `inf`.
%
%   fitScale(..., 'log', true) calculates the fit in log-space instead. See
%   the 'log' and 'logSlope' options below.
%
%   Options:
%    'mask'
%       Mask for positions in the `pred` and `meas` vectors to use when
%       finding the optimal factor `a`.
%    'stds'
%       Error bars for each measurement. This can be a single vector, or a
%       2-column matrix giving the low and high values for each point. In
%       the latter case, the function does the right thing to calculate the
%       standard deviation depending on the choice for the 'log' option
%       below.
%    'log'
%       If true, the fit instead attempts to optimize the difference
%       between `log(predScaled)` and `log(meas)`, where
%           predScaled = a*pred .
%       In this case, the errors are assumed to be distributed according to
%       a log-normal, so that the objective function that is to be
%       minimized (and that appears in the return argument `mse`) is
%           mse = mean( ((log(predscaled) - log(meas)) ./ stds) .^2 )
%    'logSlope'
%       If true, the fit uses both a slope and an intercept in log space,
%           predscaled = a(1)*pred.^a(2) .
%    'fixScale'
%       Instead of fitting, fix the scale (`a`) to the value(s) given here
%       and calculate the scaled predictions and the MSE.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% prepare to show defaults
% (doing this because some defaults depend on `meas` being provided)
if nargin == 1 && strcmp(pred, 'defaults')
    showDefaults = true;
    meas = [];
else
    showDefaults = false;
end

parser.addParameter('stds', [], @(v) isnumeric(v) && ((isvector(v) && length(v) == length(meas)) || ...
    (ismatrix(v) && size(v, 1) == length(meas) && size(v, 2) == 2)));
parser.addParameter('mask', true(size(meas)), @(b) islogical(b) && isvector(b) && length(b) == length(meas));
parser.addParameter('log', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('logSlope', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('fixScale', [], @(a) isnumeric(a) && isvector(a) && ismember(length(a), [1 2]));

% show defaults
if showDefaults
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
elseif ~isvector(params.stds)
    if ~params.log
        params.stds = diff(params.stds, [], 2);
    else
        params.stds = diff(log(params.stds), [], 2);
    end
end

% figure out which entries to care about
measMask = isfinite(meas) & isfinite(params.stds) & params.mask;
predMask = isfinite(pred) & params.mask;

if params.log
    % need the values to also be positive
    measMask = measMask & (meas > 0);
    predMask = predMask & (pred > 0);
end

% these are the common positions where both predictions and measurements
% are valid
mask = measMask & predMask;

if ~params.log
    % dividing by error bars maps to problem where error bars are constant
    predNorm = pred(mask) ./ params.stds(mask);
    measNorm = meas(mask) ./ params.stds(mask);

    % find scaling factor
    if isempty(params.fixScale)
        a = dot(predNorm, measNorm) / dot(predNorm, predNorm);
    else
        a = params.fixScale;
    end
    predScaled = a*pred;
else
    % doing the fit in log space
    x = log(pred(mask));
    y = log(meas(mask));
    logWeights = 1./params.stds(mask).^2;
    if ~params.logSlope
        % things are easy if we don't fit a slope in log space
        % (i.e., we assumed predictions are just a scaled version of the
        %  measurements, with no power-law warping)
        logA = wmean(y - x, logWeights);
        % xScaled is used for calculating the MSE
        xScaled = x + logA;
        % set the output arguments
        if isempty(params.fixScale)
            a = exp(logA);
        else
            a = params.fixScale;
        end
        predScaled = a*pred;
    else
        % we use both a slope and an intercept
        % precalculate some values we use repeatedly
        xwm = wmean(x, logWeights);
        ywm = wmean(y, logWeights);
        xywm = wmean(x.*y, logWeights);
        x2wm = wmean(x.^2, logWeights);
        
        % the scaling factor output is now a 2-component vector
        if isempty(params.fixScale)
            a(2) = (xywm - xwm*ywm) / (x2wm - xwm^2);
        else
            a = params.fixScale;
        end
        logA = wmean(y - a(2)*x, logWeights);
        
        % xScaled is used for calculating the MSE
        xScaled = a(2)*x + logA;
        
        if isempty(params.fixScale)
            a(1) = exp(logA);
        end
        
        % when applying the transformation to the predictions, we can only
        % do it for positive entries
        posMask = (pred > 0);
        predScaled = nan(size(pred));
        predScaled(posMask) = a(1)*pred(posMask).^a(2);
    end
end

% don't bother calculating the MSE unless it's actually used
if nargout > 2
    if sum(mask) > 0
        % if we don't have predictions for any of the entries where we have
        % good measurements, we're doing a bad fit
        if sum(measMask & ~predMask) > 0
            mse = inf;
        else
            if ~params.log
                mse = mean( ((predScaled(mask) - meas(mask)) ./ params.stds(mask)) .^2 );
            else
                mse = wmean((xScaled - y) .^ 2, logWeights);
            end
        end
    else
        % we don't have any common predictions and measurements
        mse = nan;
    end
end

end

function m = wmean(v, w)
% WMEAN Weighted mean.

m = mean(v .* w) / mean(w);

end