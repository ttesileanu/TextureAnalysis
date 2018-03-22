function L = solveLinearEfficientCoding(S, covNoiseIn, covNoiseOut, lag, varargin)
% solveLinearEfficientCoding Find the optimal gains matrix from a linear
% efficient-coding argument.
%   L = solveLinearEfficientCoding(S, covNoiseIn, covNoiseOut, lag)
%   calculates the optimal linear transformation that encodes a
%   multivariate Gaussian-distributed input signal with covariance matrix
%   `S` into an output signal of dimension `size(covOutNoise, 1)`, assuming
%   a constraint on the output power related to the Lagrange multiplier
%   `lag`. The signal is assumed to be corrupted by input noise with
%   covariance matrix `covNoiseIn`, and the processed signal is assumed to
%   be corrupted by output signal with covariance matrix `covNoiseOut`. If
%   `m = size(covNoiseOut, 1)` and `n = size(covNoiseIn, 1)`, the output
%   matrix `L` will have size `(m, n)`.
%
%   XXX Currently, covNoiseOut must be a multiple of the identity for the
%   result to be correct.
%
%   To obtain non-zero gains, the Lagrange multiplier `lag` should be
%   between 0 and max(1./eig(covNoiseOut)).
%
%   Options:
%    'tol'
%       Tolerance level for checking symmetry of covariance matrices. For a
%       covariance matrix A, if max(abs(flatten(A - A'))) <= tol, the matrix
%       will be considered symmetric, and replaced by (A + A')/2.
%       Otherwise, an error message will be issued.
%    'covreg'
%       To avoid problems with eigenvalues that are negative but extremely
%       tiny, we add a multiple of identity to the input noise covariance
%       matrix before taking the square root. This parameter sets the
%       multiple.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('tol', 1e-6, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('covreg', 1e-8, @(x) isscalar(x) && isnumeric(x));

if nargin == 1 && strcmp(S, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% check the data: make sure Lagrange multiplier is scalar
if ~isscalar(lag)
    error([mfilename ':lagnonsc'], 'Lagrange multiplier `lag` should be a scalar.');
end

% check the data: make sure inputs are square matrices
checkSquareShape(S, 'S');
checkSquareShape(S, 'covNoiseIn');
checkSquareShape(S, 'covNoiseOut');

% check the data: make sure dimensions agree
if ~all(size(S) == size(covNoiseIn))
    error([mfilename ':sizemismatch'], 'The sizes of the signal and input noise covariance matrices don''t agree.');
end

m = size(covNoiseOut, 1);
n = size(covNoiseIn, 1);

% XXX for now we can only handle reduced dimensionality
if m > n
    error([mfilename ':notimp'], 'Currently only cases in which the input has at least as many dimensions as the output are supported.');
end

% check the data: make sure covariance matrices are symmetric
S = checkSymmetry(S, params.tol, 'Signal');
covNoiseIn = checkSymmetry(covNoiseIn, params.tol, 'Input noise');
covNoiseOut = checkSymmetry(covNoiseOut, params.tol, 'Output noise');

% rotate away input noise by using its matrix square root
sqrtNoiseIn = sqrtm(covNoiseIn + params.covreg*eye(size(covNoiseIn)));
invSqrtNoiseIn = inv(sqrtNoiseIn);

Shat = invSqrtNoiseIn*S*invSqrtNoiseIn'; %#ok<MINV>

% move to the principal component coordinate system for the signal
[pcS, diagS] = eigsorted(Shat);
signalVars = diag(diagS);

% rotate away output noise by using its eigendecomposition
[evecNoiseOut, diagNoiseOut] = eig(covNoiseOut);
noiseOutVars = diag(diagNoiseOut);

% XXX currently I'm not sure how to do the case of non-spherically-symmetric
% output noise
if std(noiseOutVars) > params.tol
    error([mfilename ':notimp'], 'Currently only spherically-symmetric output noise is supported.');
end

% calculate gains in signal-PC coordinates
rotatedGainsSquared = zeros(m, 1);
discriminants = signalVars.^2 + (4/lag)*signalVars./noiseOutVars;
mask = discriminants >= 0;
rotatedGainsSquared(mask) = noiseOutVars(mask) .* ...
    (-(2 + signalVars(mask)) + sqrt(discriminants(mask))) ./ (2*(1 + signalVars(mask)));

rotatedGains = sqrt(max(rotatedGainsSquared, 0));

Lrotated= full(spdiags(rotatedGains, 0, m, n));

L = evecNoiseOut*Lrotated*pcS'*invSqrtNoiseIn; %#ok<MINV>

end

function A = checkSymmetry(A, tol, name)
% Check whether the matrix is symmetric with given tolerance. If it is,
% return (A + A')/2. If not, throw an error, using the given name for the
% matrix.

if ~isempty(name)
    name(1) = upper(name(1));
end

if max(abs(flatten(A - A'))) > tol
    error([mfilename ':nonsym'], [name ' covariance matrix is not symmetric.']);
else
    A = (A + A')/2;
end

end

function checkSquareShape(A, name)
% Check whether the input argument is a square matrix, and if it's not,
% issue an error.

if ~ismatrix(A)
    error([mfilename ':notmat'], [name ' should be a matrix.']);
end

if size(A, 1) ~= size(A, 2)
    error([mfilename ':notsquare'], [name ' should be square.']);
end

end