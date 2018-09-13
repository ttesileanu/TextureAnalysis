function L = solveLinearEfficientCoding(S, varargin)
% solveLinearEfficientCoding Find the optimal gains matrix from a linear
% efficient-coding argument.
%   L = solveLinearEfficientCoding(S) calculates the optimal linear
%   transformation that encodes a multivariate Gaussian-distributed input
%   signal with covariance matrix `S`, assuming isotropic input and output
%   noise, and also assuming that the output noise is dominant over *both*
%   the input noise and the signal. This is the "variance = salience"
%   regime of efficient coding. See the options for more general cases.
%
%   Options:
%    'covNoiseIn'
%       Control the covariance matrix of the input noise.
%    'covNoiseOut'
%       Control the covariance matrix of the output noise. This can also be
%       used to indicate that the dimensionality of the output should be
%       different from that of the input.
%       XXX Currently the function only works with isotropic output noise.
%       XXX For now the function can only handle reduced dimensionality.
%    'lagrange'
%       If this is provided, the function uses more general efficient
%       coding calculations that do not assume the "variance = salience"
%       regime. The value provided here is used as the Lagrange multiplier
%       related to the constraint on output power from the filter. To obtain
%       non-zero gains, the Lagrange multiplier should be between 0 and
%       `max(1./eig(covNoiseOut))`.
%    'gainTransform'
%       A transformation function to apply to gains in the principal
%       component coordinates.
%    'tol'
%       Tolerance level for checking symmetry of covariance matrices. For a
%       covariance matrix A, if max(abs(flatten(A - A'))) <= tol, the matrix
%       will be considered symmetric, and replaced by (A + A')/2. Otherwise,
%       an error message will be issued.
%    'covarianceRegularization'
%       To avoid problems with eigenvalues that are negative but extremely
%       tiny, we add a multiple of identity to the input noise covariance
%       matrix before taking the square root. This parameter sets the
%       multiple.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('covNoiseIn', [], @(m) ismatrix(m) && isscalar(m));
parser.addParameter('covNoiseOut', [], @(m) ismatrix(m) && isscalar(m));
parser.addParameter('lagrange', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
parser.addParameter('gainTransform', @(x) x, @(f) isa(f, 'function_handle'));
parser.addParameter('tol', 1e-6, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('covarianceRegularization', 1e-8, @(x) isscalar(x) && isnumeric(x));

% show defaults
if nargin == 1 && strcmp(S, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% set defaults
if isempty(params.covNoiseIn)
    params.covNoiseIn = eye(size(S));
else
    if ~isequal(size(params.covNoiseIn), size(S))
        error([mfilename ':badcovin'], 'Size of the input noise covariance matrix doesn''t match the input signal size.');
    end
end
if isempty(params.covNoiseOut)
    params.covNoiseOut = eye(size(S));
else
    % XXX for now we can only handle reduced dimensionality
    if size(params.covNoiseOut, 1) > size(params.covNoiseIn, 1)
        error([mfilename ':notimp'], 'Currently only cases in which the input has at least as many dimensions as the output are supported.');
    end
end

% check the data: make sure inputs are square matrices
checkSquareShape(S, 'S');
checkSquareShape(params.covNoiseIn, 'covNoiseIn');
checkSquareShape(params.covNoiseOut, 'covNoiseOut');

m = size(params.covNoiseOut, 1);
n = size(params.covNoiseIn, 1);

% check the data: make sure covariance matrices are symmetric
S = checkSymmetry(S, params.tol, 'Signal');
covNoiseIn = checkSymmetry(params.covNoiseIn, params.tol, 'Input noise');
covNoiseOut = checkSymmetry(params.covNoiseOut, params.tol, 'Output noise');

% rotate away input noise by using its matrix square root
sqrtNoiseIn = sqrtm(covNoiseIn + params.covarianceRegularization*eye(size(covNoiseIn)));
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

% calculate gains-squared in signal-PC coordinates
rotatedGainsSquared = zeros(m, 1);
if ~isempty(params.lagrange)
    discriminants = signalVars.^2 + (4/params.lagrange)*signalVars./noiseOutVars;
    mask = discriminants >= 0;
    rotatedGainsSquared(mask) = noiseOutVars(mask) .* ...
        (-(2 + signalVars(mask)) + sqrt(discriminants(mask))) ./ (2*(1 + signalVars(mask)));
else
    % "variance = salience" approximation
    discriminants = noiseOutVars .* signalVars;
    mask = discriminants >= 0;
    rotatedGainsSquared(mask) = sqrt(discriminants(mask));
end

% get the gains themselves, and rotate back to the original coordinates
rotatedGains = sqrt(max(rotatedGainsSquared, 0));

% apply the transform
rotatedGains = params.gainTransform(rotatedGains);

Lrotated = full(spdiags(rotatedGains, 0, m, n));
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

if ~isempty(name)
    name(1) = upper(name(1));
end

if ~ismatrix(A)
    error([mfilename ':notmat'], [name ' should be a matrix.']);
end

if size(A, 1) ~= size(A, 2)
    error([mfilename ':notsquare'], [name ' should be square.']);
end

end
