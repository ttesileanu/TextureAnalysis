function ellipse(x0, y0, a, b, theta, varargin)
% ELLIPSE Draw an ellipse.
%   ELLIPSE(x0, y0, a, b, theta) draws an ellipse centered at x0, y0, with
%   semiaxes a, b, such that the 'a' semiaxis makes an angle theta with the
%   x-axis (in the clockwise direction).
%
%   ELLIPSE(x0, y0, M), where M is a 2x2 matrix, uses the eigenvalues and
%   eigenvectors of the symmetric, positive-definite M matrix to draw the
%   ellipse. This is equivalent to drawing the locus of points r where
%   r'*M*r = 1.
%
%   Additional arguments are directly passed to Matlab's PLOT, except for:
%    'npoints'
%       Set the number of points used for drawing the ellipse.
%
%   See defaults by running ELLIPSE('defaults').
%
%   See also: PLOT.


% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('npoints', 48, @(x) isscalar(x) && isnumeric(x) && x >= 1);

% handle the 'defaults' option
if nargin == 1 && strcmp(x0, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% handle the two different forms of the command
if ~isscalar(a)
    varargin = [{b, theta} varargin];
    M = a;
    [V, D] = eig(M);
    D = diag(D);
    a = sqrt(D(1));
    b = sqrt(D(2));
    
    e1 = V(:, 1);
    e2 = V(:, 2);
else
    e1 = [cos(theta) ; -sin(theta)];
    e2 = [sin(theta) ; cos(theta)];
end

% separate special arguments from PLOT arguments
special_args0 = false(size(varargin));
special_args0(ischar(varargin)) = ismember(varargin(ischar(varargin)), {'npoints'});
if any(special_args0)
    if special_args0(end)
        error([mfilename ':badspargs'], 'Incomplete argument.');
    end
    special_idxs = find(special_args0);
    special_args0(special_idxs+1) = true;
    special_args = varargin(special_args0);
    varargin = varargin(~special_args0);
else
    special_args = {};
end

% parse special arguments
parser.parse(special_args{:});
params = parser.Results;

% draw the ellipse
phi = linspace(0, 2*pi, params.npoints);
c = cos(phi);
s = sin(phi);

r = a*e1*c + b*e2*s;

plot(x0 + r(1, :), y0 + r(2, :), varargin{:});

end