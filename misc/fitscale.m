function [a, predscaled, sse] = fitscale(pred, meas, varargin)
% FITSCALE Find the best scaling factor to match a vector of predictions to
% data.
%   a = FITSCALE(pred, meas) finds the factor `a` such that
%   norm(a*pred - meas) is minimized. Entries whether either predictions or
%   measurements are infinite or NaN are ignored.
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
%   [..., sse] = FITSCALE(...) also returns the summed squared error
%   obtained with the optimal scaling coefficient. Specifically,
%       sse = 1/2*sum( ((a*pred - meas) ./ stds) .^2 )
%   If error bars aren't given, `stds` is assumed to be identically 1.

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
mask = isfinite(pred) & isfinite(meas) & isfinite(params.stds) & ~params.exclude;

% dividing by error bars maps to problem where error bars are constant
prednorm = pred(mask) ./ params.stds(mask);
measnorm = meas(mask) ./ params.stds(mask);

% find scaling factor
a = dot(prednorm, measnorm) / dot(prednorm, prednorm);
predscaled = a*pred;

if nargout > 2
    if sum(mask) > 0
        sse = 1/2*sum( ((predscaled(mask) - meas(mask)) ./ params.stds(mask)) .^2 );
    else
        sse = nan;
    end
end

end