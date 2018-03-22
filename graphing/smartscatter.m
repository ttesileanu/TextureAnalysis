function [handles, colors] = smartscatter(x, y, varargin)
% SMARTSCATTER Make a scatter plot with smart subsampling for large
% datasets.
%   SMARTSCATTER(x, y) makes a scatter plot similar to Matlab's scatter
%   command. Unlike that function, SMARTSCATTER automatically flattens the
%   x and y data before plotting, and plots a subset of the data in cases
%   in which the number of points is too large. There are also easy-to-use
%   commands to control transparency (see below).
%
%   h = SMARTSCATTER(...) return handles for the scatter objects.
%
%   Key-value options that are not recognized by SMARTSCATTER are directly
%   passed to the SCATTER command.
%
%   SMARTSCATTER specific options:
%    'alpha': double
%       Transparency level of markers (0 fully transparent, 1 fully
%       opaque).
%    'density': bool
%       If true, the color of the dots is chosen based on the local density
%       of points.
%    'kernelscale':
%       Scaling to use for the kernel function. This can be a pair of
%       numbers, giving a different scale in each direction. This
%       effectively sets the range over which to blur the distribution to
%       estimate colors. By default this is chosen to be a fraction of the
%       entire data range.
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
% keep unmatched key/value pairs to pass to scatter
parser.KeepUnmatched = true;

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
    rng_x = max(x) - min(x);
    rng_y = max(y) - min(y);
    % make sure range isn't zero
    rng_x = rng_x + eps(rng_x);
    rng_y = rng_y + eps(rng_y);
    params.kernelscale = [rng_x/24, rng_y/24];
elseif length(params.kernelscale) == 1
    params.kernelscale = [params.kernelscale, params.kernelscale];
end

% get density values
% these can be used either for the density plot or for subsampling
x_scaled = x/params.kernelscale(1);
y_scaled = y/params.kernelscale(2);
% figure out the range of the data
x_min = min(x_scaled);
x_max = max(x_scaled);
y_min = min(y_scaled);
y_max = max(y_scaled);
% add some edges
edge_x = (x_max - x_min)*0.1;
edge_y = (y_max - y_min)*0.1;
x_min = x_min - edge_x;
y_min = y_min - edge_y;
% figure out locations of bin centers
stepx = (x_max - x_min)/params.densitybins;
stepy = (y_max - y_min)/params.densitybins;
locs{1} = (x_min + stepx/2):stepx:(x_max - stepx/2);
locs{2} = (y_min + stepy/2):stepy:(y_max - stepy/2);
% find bin assignments for each point
bin_x = min(max(1 + round((x_scaled - locs{1}(1))/stepx), 1), length(locs{1}));
bin_y = min(max(1 + round((y_scaled - locs{2}(1))/stepy), 1), length(locs{2}));
% generate density matrix
density_matrix0 = accumarray([bin_x, bin_y], ones(length(bin_x), 1));

if ~isempty(params.color)
    params.density = false;
end

% convolve densities with a kernel if we want to use the density for
% coloring points
if params.density    
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
    
%     % subsample only in high-density regions, to avoid missing outliers
%     max_avg_density = params.maxpoints/numel(density_matrix0);
%     high_density_mask = (density_matrix0 >= max_avg_density);
%     
%     % number of points in low-density regions -- these will always be kept
%     % (unless by an odd combination of parameters they are more than
%     %  maxpoints/2)
%     n_low_density = sum(flatten(density_matrix0(~high_density_mask)));
%     if n_low_density > params.maxpoints/2
%     if true
%         % not much that we can do: just subsample randomly
%         idxs = randperm(n, params.maxpoints);
%     else
%         % convert bin_x, bin_y to linear indices in the density matrix
%         linear_idx = sub2ind(size(density_matrix0), bin_x, bin_y);
%         % make a mask of points to keep
%         keep_mask = false(size(linear_idx));
%         n_keep_per_bin = round(max_avg_density);
%         for crt_idx = 1:numel(density_matrix0)
%             crt_mask = (linear_idx == crt_idx);
%             if high_density_mask(crt_idx)
%                 % keep first n_keep_per_bin
%                 crt_idxs_all = find(crt_mask);
%                 crt_idxs_to_keep = crt_idxs_all(randperm(length(crt_idxs_all), n_keep_per_bin));
%                 keep_mask(crt_idxs_to_keep) = true;
%             else
%                 % keep_all
%                 keep_mask(crt_mask) = true;
%             end
%         end
%         idxs = keep_mask;
%     end
    
    x = x(idxs);
    y = y(idxs);
    
    if ~isempty(params.size) && ~isscalar(params.size)
        params.size = params.size(idxs);
    end
    if ~isempty(params.color)
        if isvector(params.color) && numel(params.color) == n
            params.color = params.color(idxs);
        elseif size(params.color, 2) == 3 && size(params.color, 1) == n
            params.color = params.color(idxs, :);
        end
    end
    
    rng(rng_prev);
    
%     n = length(x);
end

if params.density
    x_scaled = x/params.kernelscale(1);
    y_scaled = y/params.kernelscale(2);
    density = interp2(locs{2}, locs{1}, density_matrix, y_scaled, x_scaled);
    params.color = params.densityfunction(density);
    % fix dots that are on the border
    params.color(~isfinite(params.color)) = min(params.color(isfinite(params.color)));
end

% figure out scatter options
if ~isempty(params.axes)
    opts = {params.axes};
else
    opts = {};
end
if ~isempty(params.size)
    opts = [opts {params.size}];
elseif ~isempty(params.color)
    opts = [opts {[]}];
end
if ~isempty(params.color)
    opts = [opts {params.color}];
end
if ~isempty(params.marker)
    opts = [opts {params.marker}];
end
if params.filled
    opts = [opts {'filled'}];
end
if params.alpha < 1
    opts = [opts {'markerfacealpha', params.alpha, 'markeredgealpha', params.alpha}];
end
opts = [opts structToCell(parser.Unmatched)];

% plot!
handles.hscatter = scatter(x, y, opts{:});

colors = params.color;

end