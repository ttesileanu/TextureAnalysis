function [intersects, mapped, details] = mapIntersectEllipsoid(interpolation, trafo, radius2, varargin)
% mapIntersectEllipsoid Intersect the rays of an interpolated mapping with
% an ellipsoid.
%   [intersects, mapped] = mapIntersectEllipsoid(interpolation, trafo, radius2)
%   goes through all the functions mapping scalar keys to vector values from
%   the cell array of structures `interpolation` (as returned by
%   mapInterpolate) and finds the point where each of these mappings
%   intersect the ellipsoid given by
%
%       v(:)'*trafo*v(:) = radius2 .
%
%   The resulting keys are returned in the `intersects` vector, while the
%   corresponding positions mapped through the functions from `interpolation`
%   are returned in the `mapped` cell array.
%
%   [..., details] = mapIntersectEllipsoid(...) also returns a cell array
%   containing the output from FSOLVE for each direction. This is given in
%   a structure containing the `fval`, `exitflag`, `output`, and `jacob`
%   output arguments of FSOLVE.
%
%   Options:
%    'guess'
%       Guess value(s) to use for finding the intersection point. This can
%       be a scalar, in which case the same guess value is used for every
%       direction, or a vector giving the guess value for each direction.
%    'fsolve_opts'
%       Extra options to pass to FSOLVE. The same options are used for all
%       directions.
%    'force_positive'
%       If true, the intersects are forced to be positive.
%
%   See also: FSOLVE.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% generate default fsolve options
default_opts = optimoptions('fsolve', 'display', 'none');

% add possible parameters with their defaults
parser.addParameter('guess', 0.1, @(v) isnumeric(v) && isvector(v));
parser.addParameter('fsolve_opts', {default_opts}, @(c) iscell(c));
parser.addParameter('force_positive', true, @(b) islogical(b) && isscalar(b));

% handle showing of defaults.
if nargin == 1 && strcmp(interpolation, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse optional arguments
parser.parse(varargin{:});
params = parser.Results;

% handle scalar vs. vector guess values
if isscalar(params.guess)
    params.guess = repmat(params.guess, size(interpolation));
end

% norm using the trafo matrix
trafo_norm = @(v) v(:)'*trafo*v(:);

% initialize outputs
intersects = zeros(size(interpolation));
mapped = cell(size(interpolation));
details = cell(size(interpolation));

% go through each direction, find intersection
for i = 1:length(interpolation)
    crt_fct = interpolation{i}.function;
    crt_guess = params.guess(i);
    
    if ~params.force_positive
        [crt_x, crt_fval, crt_exitflag, crt_output, crt_jacob] = fsolve(...
            @(x) trafo_norm(crt_fct(x)) - radius2, crt_guess, ...
            params.fsolve_opts{:});
    else
        [crt_sqrtx, crt_fval, crt_exitflag, crt_output, crt_jacob] = fsolve(...
            @(sqrt_x) trafo_norm(crt_fct(sqrt_x^2)) - radius2, sqrt(crt_guess), ...
            params.fsolve_opts{:});
        crt_x = crt_sqrtx^2;
    end

    intersects(i) = crt_x;
    mapped{i} = crt_fct(crt_x);
    
    crt_struct = struct;
    crt_struct.fval = crt_fval;
    crt_struct.exitflag = crt_exitflag;
    crt_struct.output = crt_output;
    crt_struct.jacob = crt_jacob;
    details{i} = crt_struct;
end

end