function showTernaryTexturesSingle(group, varargin)
% showTernaryTexturesSingle Show example textures in single texture plane.
%   showTernaryTexturesSingle(group) displays example texture patches at
%   the corners of the texture simplex, drawing samples from the given
%   texture `group`.
%
%   Options:
%    'patchSize'
%       Size of the patches in number of checks for each dimension.
%    'patchExtent'
%       Size of the patches in data coordinates.
%    'offset'
%       Amount by which the center of the patches is displaced with respect
%       to the vertex that they correspond to, as a fraction of the axis
%       length.
%    'bimageOpts'
%       Options to be passed to `bimage` for drawing the patches.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('patchSize', 12, @(n) isscalar(n) && isnumeric(n) && n > 0);
parser.addParameter('patchExtent', 0.7, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('offset', 0.5, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('bimageOpts', {'borderwidth', 0.5}, @(c) iscell(c));

% show defaults if requested
if nargin == 1 && strcmp(group, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% set hold, but maintain old state
wasHold = ishold;
hold on;

% generate and draw the patches
patches = cell(1, 3);
for i = 1:3
    % set up the direction
    direction = zeros(1, 3);
    direction(i) = 1;
    
    % set up the generator
    generator = PatchAxisGenerator(group, direction, params.patchSize);
    generator.locations = 1;
    
    % generate
    generator.next;
    patches{i} = generator.samples;
    
    % find the location of the vertex in the plane
    vertex2 = ternary3to2(direction);
    
    % calculate the location of the upper-left corner of the patch
    
    % need slightly larger offset in diagonal directions
    fudge = (all(abs(vertex2) > 0))*0.2;
    center = vertex2*(1 + params.offset*(1 + fudge));
    corner = center - 0.5*[params.patchExtent params.patchExtent];
    
    % draw the patch
    bimage([corner(1) corner(1) + params.patchExtent], ...
           [corner(2) + params.patchExtent corner(2)], ...
           patches{i}, 'CDataMapping', 'scaled');
end
set(gca, 'clim', [0 1]);
colormap('gray');

if ~wasHold
    hold off;
end

end
