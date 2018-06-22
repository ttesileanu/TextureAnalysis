function [handles, colors] = smartscatter(x, y, varargin)
% SMARTSCATTER Make a scatter plot with smart subsampling for large
% datasets and possibility for density-based coloring.
%   SMARTSCATTER(x, y) makes a scatter plot similar to Matlab's scatter
%   command. Unlike that function, SMARTSCATTER automatically flattens the
%   x and y data before plotting, and plots a subset of the data in cases
%   in which the number of points is too large. The function also
%   automatically ignores entries for which either of the input vectors
%   is infinite or NaN. There are also easy-to-use commands to control
%   transparency (see below) and coloring of the dots based on their local
%   density.
%
%   h = SMARTSCATTER(...) return handles for the scatter objects.
%
%   [h, colors] = SMARTSCATTER(...) returns the final color vector that was
%   used. This is useful when 'density' is true.
%
%   Options:
%    'alpha': double
%       Transparency level of markers (0 fully transparent, 1 fully
%       opaque).
%    'density': bool
%       If true, the color of the dots is chosen based on the local density
%       of points.
%    'kernelscale':
%       Scaling to use for the kernel function when using density coloring.
%       This can be a pair of numbers, giving a different scale in each
%       direction. This effectively sets the range over which to blur the
%       distribution to estimate colors. By default this is chosen to be a
%       fraction of the entire data range.
%    'densityrange'
%       How far to look for neighboring points when counting density, in
%       units of `kernelscale`. Too large a range can make density
%       estimation very slow.
%    'densitybins'
%       Number of bins to use in each direction for density estimation.
%    'densityfunction'
%       A function to apply to the density values before setting them as
%       scatter colors.
%    'maxpoints' <n>
%       Maximum number of points to display. If the number of elements in
%       the input vectors is larger than this, a uniform sampling will be
%       done so that only n points are plotted.
%    'size': vector
%       Vector of marker sizes. If this is not a vector, it is flattened
%       before drawing.
%    'color': vector, matrix
%       Vector or matrix of marker colors. This can be a vector of scalars
%       that will get mapped to colors using the current color map. Or it
%       can be an N x 3 matrix of RGB colors.
%    'axes':
%       Axis where to draw the scatterplot, if different from gca.
%    'filled':
%       Whether to use filled markers or not.
%    'marker':
%       Marker to use.
%
%   SMARTSCATTER('defaults') shows the defaults for these options.
%
%   See also: SCATTER.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('alpha', 0.6, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('density', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('maxpoints', 10000, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('size', [], @(v) isnumeric(v) && (isscalar(v) || numel(v) == numel(x)));
parser.addParameter('color', [], @(v) isnumeric(v) && ((isvector(v) && ...
    numel(v) == numel(x)) || (ismatrix(v) && size(v, 2) == 3 && ...
    (size(v, 1) == numel(x) || size(v, 1) == 1))));
parser.addParameter('axes', []);
parser.addParameter('filled', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('marker', [], @(s) ischar(s) && isscalar(s));
parser.addParameter('kernelscale', [], @(s) isnumeric(s) && isvector(s) && ...
    length(s) <= 2);
parser.addParameter('densityrange', 1, @(s) isnumeric(s) && isscalar(s));
parser.addParameter('densityfunction', @(v) v, @(f) isa(f, 'function_handle'));
parser.addParameter('densitybins', 30, @(n) isnumeric(n) && isscalar(n) && n > 1);

if nargin == 1 && strcmp(x, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% flatten inputs
x = x(:);
y = y(:);

% check that the sizes match
n = length(x);
if length(y) ~= n
    error([mfilename ':badsz'], 'Numbers of elements of x and y must match.');
end

% get rid of NaNs / infinities
mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);

% handle defaults
if isempty(params.kernelscale)
    % automatic kernel scale, as a fraction of the data range
    rng_x = max(x) - min(x);
    rng_y = max(y) - min(y);
    % make sure range isn't zero
    rng_x = rng_x + eps(rng_x);
    rng_y = rng_y + eps(rng_y);
    params.kernelscale = [rng_x/24, rng_y/24];
elseif isscalar(params.kernelscale)
    params.kernelscale = [params.kernelscale, params.kernelscale];
end

% if an explicit color is provided, disable density plot
if ~isempty(params.color)
    params.density = false;
end

% for making density plot, calculate densities in a grid with cell size
% given by the kernelsize, then smooth by convolving with a kernel
if params.density
    % get raw density values
    x_scaled = x/params.kernelscale(1);
    y_scaled = y/params.kernelscale(2);
    % figure out the range of the data
    x_min = min(x_scaled);
    x_max = max(x_scaled);
    y_min = min(y_scaled);
    y_max = max(y_scaled);
    % choose bin sizes
    stepx = (x_max - x_min)/(params.densitybins - 1);
    stepy = (y_max - y_min)/(params.densitybins - 1);
    % figure out locations of bin centers
    locs{1} = x_min:stepx:x_max;
    locs{2} = y_min:stepy:y_max;
    % find bin assignments for each point
    bin_x = min(max(1 + round((x_scaled - locs{1}(1))/stepx), 1), length(locs{1}));
    bin_y = min(max(1 + round((y_scaled - locs{2}(1))/stepy), 1), length(locs{2}));
    % generate density matrix
    density_matrix0 = accumarray([bin_x, bin_y], ones(length(bin_x), 1));

    [kx, ky] = meshgrid(-3:stepx:3, -3:stepy:3);
    kernel_matrix = exp(-(kx.^2 + ky.^2)/2)/(2*pi);
    density_matrix = conv2(density_matrix0, kernel_matrix, 'same');
end

% subsample if necessary
if n > params.maxpoints
    % when sampling randomly, do it reproducibly
    rng_prev = rng;
    rng('default');
    
    % first thought I'd only subsample high-density areas, but that doesn't
    % work well: either we cap the density to a fixed value, but then we
    % lose the sense of density from the scatter plot; or we subsample
    % randomly only in high-density areas, leaving low-density areas as
    % they are (to see outliers), but then there is a gap in the transition
    % between low density and high density
    idxs = randperm(n, params.maxpoints);
    
    x = x(idxs);
    y = y(idxs);
    
    % subsample the sizes if needee
    if ~isempty(params.size) && ~isscalar(params.size)
        params.size = params.size(idxs);
    end
    % subsample the colors if needed
    if ~isempty(params.color)
        if isvector(params.color) && numel(params.color) == n
            params.color = params.color(idxs);
        elseif size(params.color, 2) == 3 && size(params.color, 1) == n
            params.color = params.color(idxs, :);
        end
    end
    
    % go back to the previous random number generator state
    rng(rng_prev);
    
%     n = length(x); % <-- not using n anymore
end

% assign colors based on point density, if necessary
if params.density
    % might have subsampled x&y, so just recalculate scaled versions
    x_scaled = x/params.kernelscale(1);
    y_scaled = y/params.kernelscale(2);
    % use bilinear interpolation to find density at each point
    density = interp2(locs{2}, locs{1}, density_matrix, y_scaled, x_scaled);
    % apply a function to the density -- useful to compress densities with
    % high dynamic range
    params.color = params.densityfunction(density);
    % fix any dots that might have been outside the interpolation range
    params.color(~isfinite(params.color)) = min(params.color(isfinite(params.color)));
end

% figure out scatter options

% use specific axes?
if ~isempty(params.axes)
    opts = {params.axes};
else
    opts = {};
end

% the data to scatter
opts = [opts {x y}];

% size information?
if ~isempty(params.size)
    opts = [opts {params.size}];
elseif ~isempty(params.color)
    opts = [opts {[]}];
end

% color information?
if ~isempty(params.color)
    opts = [opts {params.color}];
end

% set the marker
if ~isempty(params.marker)
    opts = [opts {params.marker}];
end

% filled or not?
if params.filled
    opts = [opts {'filled'}];
end

% transparency
if params.alpha < 1
    opts = [opts {'markerfacealpha', params.alpha, 'markeredgealpha', params.alpha}];
end

% plot!
handles.hscatter = scatter(opts{:});

colors = params.color;

end