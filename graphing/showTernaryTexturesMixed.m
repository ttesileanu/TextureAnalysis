function showTernaryTexturesMixed(group, varargin)
% showTernaryTexturesMixed Show example textures in mixed texture plane.
%   showTernaryTexturesMixed(group) displays example texture patches at
%   several positions around a mixed texture plane, identified by the
%   `group` argument. This should have the format
%       'group1[dir1];group2[dir2]'
%   where `group1` and `group2` should be valid group names, and `dir1` and
%   `dir2` should be integers from 0 to 2, identifying the direction in
%   each plane.
%
%   Options:
%    'patchSize'
%       Size of the patches in number of checks for each dimension.
%    'patchExtent'
%       Size of the patches in data coordinates.
%    'saturation'
%       Distance from the origin of texture space where to sample the
%       texture patches.
%    'saturationScale'
%       This can be
%        'absolute': the saturation is expressed in absolute units
%        'relative': the saturation is expressed in terms of the maximum
%                    saturation achievable by the texture generation
%                    process in each direction
%    'fixDistance'
%       If non-empty, all patches are drawn this distance away from center.
%       Otherwise, the distance is given by the texture coordinates.
%    'n'
%       Number of example patches to draw. These will be distributed
%       uniformly around the circle of radius given by the saturation,
%       starting on the horizontal axis.
%    'bimageOpts'
%       Options to be passed to `bimage` for drawing the patches.
%    'drawRadii'
%       Draw radial lines from center to each of the patches.
%    'radiiOpts'
%       Options to pass to `plot` when drawing the radial lines.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('patchSize', 16, @(n) isscalar(n) && isnumeric(n) && n > 0);
parser.addParameter('patchExtent', 0.2, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('saturation', 0.33, @(x) isscalar(x) && isnumeric(x) && x >= 0 && x <= 1);
parser.addParameter('saturationScale', 'absolute', @(s) ismember(s, {'absolute', 'relative'}));
parser.addParameter('fixDistance', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x >= 0));
parser.addParameter('n', 8, @(n) isscalar(n) && isnumeric(n) && n >= 1);
parser.addParameter('bimageOpts', {'borderwidth', 0.5}, @(c) iscell(c));
parser.addParameter('drawRadii', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('radiiOpts', {':', 'color', [0.5, 0.5, 0.5], 'linewidth', 0.5}, ...
    @(c) iscell(c));

% show defaults if requested
if nargin == 1 && strcmp(group, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% parse the group names
groups = cellfun(@strtrim, strsplit(group, ';'), 'uniform', false);
if length(groups) ~= 2
    error([mfilename ':badgrp'], 'This function only works for pairs of ternary groups.');
end
% parse directions in each group
directions = zeros(1, 2);
for i = 1:2
    crtIdx = find(groups{i} == '[', 1);
    if isempty(crtIdx) || crtIdx+1 > length(groups{i})
        error([mfilename ':nodir'], 'Missing direction from group %s.', groups{i});
    end
    directions(i) = str2double(groups{i}(crtIdx+1));
    groups{i} = strtrim(groups{i}(1:crtIdx-1));
end

% set hold, but maintain old state
wasHold = ishold;
hold on;

% generate and draw the patches
patches = cell(1, params.n);
for i = 1:params.n
    % set up the direction in 2d
    crtAngle = (i-1)/params.n*2*pi;
    vertex2norm = [cos(crtAngle) sin(crtAngle)];
    
    % convert to 6d direction
    direction = reshape(ternarymix2to6(vertex2norm, directions), 3, 2)';
    
    % figure out the maximum location, if necessary
    if strcmp(params.saturationScale, 'relative')
        generator = PatchAxisGenerator(groups, direction, params.patchSize);
        generator.nLocations = 1;
        generator.next;
        
        maxLocation = generator.getMaxLocation();
    else
        maxLocation = 1;
    end
    
    % set up the generator
    generator = PatchAxisGenerator(groups, direction, params.patchSize);
    generator.locations = params.saturation*maxLocation;
    
    % generate
    generator.next;
    patches{i} = generator.samples;
    
    if isempty(patches{i})
        warning([mfilename ':failgen'], ['Texture patch generation failed for ' ...
            'patch %d. Try lowering the saturation parameter.'], i);
    end
    
    % calculate the location of the upper-left corner of the patch
    if isempty(params.fixDistance)
        vertex2 = params.saturation*maxLocation*vertex2norm;
    else
        vertex2 = params.fixDistance*vertex2norm;
    end
    corner = vertex2 - 0.5*[params.patchExtent params.patchExtent];
    
    % draw radial line, if asked to
    if params.drawRadii
        plot([0 vertex2(1)], [0 vertex2(2)], params.radiiOpts{:});
    end
    
    % draw the patch
    bimage([corner(1) corner(1) + params.patchExtent], ...
           [corner(2) + params.patchExtent corner(2)], ...
           patches{i}, 'CDataMapping', 'scaled', params.bimageOpts{:});
end
set(gca, 'clim', [0 1]);
colormap('gray');

if ~wasHold
    hold off;
end

end
