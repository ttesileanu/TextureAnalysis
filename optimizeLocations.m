function locs = optimizeLocations(centers, sigmas, varargin)
% optimizeLocations Find optimal positioning of a set of labels given a set
% of ellipses.
%   locs = optimizeLocations(centers, sigma) finds the best locations for a
%   set of labels by minimizing the distances from the `centers` (a matrix
%   of size [N, 2]), and maximizing inter-label distances. The distances
%   from the centers are normalized using the standard deviations from
%   `sigmas` (a matrix of size [N, 2]), allowing for elliptically-shaped
%   centers of interest.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('labelSize', 100, @(x) isscalar(x) && isnumeric(x) && isreal(x));
parser.addParameter('repulsion', 0.2, @(x) isscalar(x) && isnumeric(x) && isreal(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% need to ignore NaNs
locs = nan(size(centers));
mask = ~any(isnan([centers sigmas]), 2);

centersValid = centers(mask, :);
sigmasValid = sigmas(mask, :);

N = sum(mask);
rN = params.repulsion*N;
varNorm = sigmasValid.^2/params.labelSize^2;

sums = sum((centersValid + params.repulsion*(N-1)/params.labelSize) ./ ...
    (1 + rN*varNorm), 1) ./ ...
    (1 - params.repulsion*sum(varNorm./(1 + rN*varNorm), 1));
locsValid = zeros(size(centersValid));
locsValid(:, 1) = ...
    (centersValid(:, 1) + params.repulsion*sums(1)*varNorm(:, 1) + params.repulsion*(N-1)/params.labelSize) ./ ...
    (1 + rN*varNorm(:, 1));
locsValid(:, 2) = ...
    (centersValid(:, 2) + params.repulsion*sums(2)*varNorm(:, 2) + params.repulsion*(N-1)/params.labelSize) ./ ...
    (1 + rN*varNorm(:, 2));

locs(mask, :) = locsValid;

% % center of mass of centers of interest
% com = mean(centers(mask, :), 1);
% 
% 
% locs(mask, :) = (centers(mask, :) - rN*sigmas(mask, :).^2/...
%     (params.labelSize^2 + rN).*(com + rN*sum(varNorm))) ./ (1 - rN*varNorm);

end