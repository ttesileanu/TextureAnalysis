function measurements = inventMeasurements(groups, directions, varargin)
% inventMeasurements Generate artificial psychophysics thresholds for
% texture sensitivity.
%   measurements = inventMeasurements(groups, directions) generates
%   thresholds in the given texture directions from a random distribution
%   controlled by the options.
%
%   This currently only works with ternary textures.
%
%   Options:
%    'method'
%       This can be
%        'independent'
%           Thresholds for each axis are drawn from log-normal
%           distributions with parameters given by the 'axisParams' option.
%           There are no correlations between axes.
%        'correlated'
%           XXX The logarithms of the thresholds are drawn from a multivariate
%           normal in which there are correlations between different axes
%           as given by the 'axisCorr' option. The distributions for each
%           axis are still given by the 'axisParams' option.
%    'axisParams'
%       Parameters defining the distribution from which thresholds are
%       drawn for each axis. These are the mean and standard deviation of
%       the underlying normal distribution.
%    'axisCorr'
%       Constant correlation between different axes. When this is set to 0,
%       the 'correlated' and 'independent' methods are the same.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('method', 'independent', @(s) ismember(s, {'independent', 'correlated'}));
parser.addParameter('axisParams', [0 1], @(v) isvector(v) && length(v) == 2 && isnumeric(v));
parser.addParameter('axisCorr', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('originalCov', [], @(m) ismatrix(m) && isnumeric(m));

% show defaults if requested
if nargin == 1 && strcmp(directions, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% find number of gray levels
G = length(directions{1}) / (1 + sum(groups{1} == ';'));

if G ~= 3
    erorr([mfilename ':notimp'], 'This function only works with ternary textures.');
end

% find number of texture groups
nGroups = G*(G^2 + G - 1);

% find number of coordinates, assuming a redundant parametrization in which
% each group is defined by a G-component probability distribution
n = G*nGroups;

% generate a covariance matrix from which we'll calculate gains
switch params.method
    case 'independent'
        if isempty(params.originalCov)
            % the thresholds will end up inversely proportional to the standard
            % deviations, so we set the variances equal to minus their squares
            diagElements = lognrnd(-2*params.axisParams(1), 2*params.axisParams(2), n, 1);
            
        else
            diagElements = diag(params.originalCov);
            diagElements = diagElements(randperm(length(diagElements)));
        end
        covMat = diag(diagElements);
    case 'correlated'
        error([mfilename ':notimp'], 'Method ''correlated'' not yet implemented.');
    otherwise
        error([mfilename ':badmeth'], 'Unknown method.');
end

% project the covariance matrix to enforce the constraints that the G
% probabilities for each group add up to 1
constraints = zeros(n, G);
% assuming the groups are sequential: first G components are group 1, next
% G are group 2, etc.
for i = 1:nGroups
    constraints((i-1)*G + 1:i*G, i) = 1/sqrt(G);
end

% the constraint vectors are orthonoral, so building a projector onto the
% constraint space is trivial
Pconstraint = constraints*constraints';

% identity minus that projector then gives the projector onto the
% orthogonal subspace to the constraint space
P = eye(n) - Pconstraint;

% project the covariance matrix
% using the fact that P = P'
covMat = P*covMat*P;

% extend PP directions to the full dimensionality of texture space
directionsExt = cellfun(@(v) v - 1/3, ternaryextdir(groups, directions), ...
    'uniform', false);

% get gains
unscaledGain = solveLinearEfficientCoding(covMat, 'gainTransform', @(x) x.^2);

% calculate thresholds from the gains
thresholds = gainsToThresholds(unscaledGain, directionsExt);

% assemble the measurements structure
measurements.groups = groups;
measurements.directions = directions;
measurements.thresholds = thresholds;

end
