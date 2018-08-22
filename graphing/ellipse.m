function status = ellipse(x0, y0, a, b, theta, varargin)
% ELLIPSE Draw an ellipse.
%   ELLIPSE(x0, y0, a, b, theta) draws an ellipse centered at `x0`, `y0`,
%   with semiaxes `a`, `b`, such that the `a` semiaxis makes an angle theta
%   with the x-axis (in the clockwise direction).
%
%   ELLIPSE(x0, y0, M), where `M` is a 2x2 matrix, uses the eigenvalues and
%   eigenvectors of the symmetric, positive-definite `M` matrix to draw the
%   ellipse. This is equivalent to drawing the locus of points `r` where
%   `r'*M*r = 1`.
%
%   Additional arguments are directly passed to Matlab's PLOT, except for:
%    'npoints'
%       Set the number of points used for drawing the ellipse.
%
%   See also: PLOT.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('nPoints', 48, @(x) isscalar(x) && isnumeric(x) && x >= 1);

% handle the 'defaults' option
if nargin == 1 && strcmp(x0, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% handle the two different forms of the command
% either way, we generate the directions for the two semiaxes, e1 and e2
if ~isscalar(a)
    % here we also need to calculate the lengths of the semiaxes
    varargin = [{b, theta} varargin];
    M = a;
    [V, D] = eig(M);
    D = diag(D);
    a = 1/sqrt(D(1));
    b = 1/sqrt(D(2));
    
    e1 = V(:, 1);
    e2 = V(:, 2);
else
    e1 = [cos(theta) ; -sin(theta)];
    e2 = [sin(theta) ; cos(theta)];
end

% separate special arguments from PLOT arguments
specialArgs0 = false(size(varargin));
specialArgs0(ischar(varargin)) = ismember(varargin(ischar(varargin)), {'nPoints'});
if any(specialArgs0)
    if specialArgs0(end)
        error([mfilename ':badspargs'], 'Incomplete argument.');
    end
    specialIdxs = find(specialArgs0);
    specialArgs0(specialIdxs+1) = true;
    specialArgs = varargin(specialArgs0);
    varargin = varargin(~specialArgs0);
else
    specialArgs = {};
end

% parse special arguments
parser.parse(specialArgs{:});
params = parser.Results;

% draw the ellipse
phi = linspace(0, 2*pi, params.nPoints);
c = cos(phi);
s = sin(phi);

r = a*e1*c + b*e2*s;

x = x0 + r(1, :);
y = y0 + r(2, :);
if isreal(x) && isreal(y)
    plot(x, y, varargin{:});
    status = true;
else
    % an ellipse couldn't be generated (presumably because the covariance
    % matrix provided was not positive-definite)
    status = false;
end

end
