function handles = smartscatter(x, y, varargin)
% SMARTSCATTER Make a scatter plot with label and support for non-trivial
% matching between x and y vectors.
%   SMARTSCATTER(x, y) makes a scatter plot similar to Matlab's scatter
%   command. Unlike that function, SMARTSCATTER automatically flattens the
%   x and y data before plotting, and plots a subset of the data in cases
%   in which the number of points is too large. A number of options allows
%   drawing labels next to each of the points, using a non-trivial mapping
%   between the components of x and those of y (which would allow these to
%   have different lengths), altering the size, shape, and color of the
%   points on an individual basis, and changing the framing of the plot.
%
%   Do note that when per-point options are chosen, the drawing can become
%   rather slow, as many plot commands might be issued.
%
%   handles = SMARTSCATTER(...) returns a structure of handles which can be
%   useful for modifying the appearance of the plot.
%
%   Options:
%     'alpha' <x>
%       If different from 1, use transparency.  Because of Matlab
%       limitations, this needs to be drawn differently than the other
%       options, using patches. Because of this, only circular shapes will
%       be supported, and the sizes will behave differently (especially
%       under rescaling of the figure).
%     'colors' <c/v/s>
%       This can be provided in various ways. It can be single character
%       identifying the color of the points in Matlab's format. It can be
%       an array of characters, one for each point. Alternatively, it can
%       be a column vector giving the RGB components of the color for all
%       the points, or an Nx3 matrix giving the RGB colors for every point.
%       (default: 'r')
%     'density' <b>
%       When this is true, a density plot is drawn instead of a scatter
%       plot. 'Colors', 'maxpoints', 'shapes', and 'sizes' are ignored, and
%       instead an image with resolution given by 'resolution' is drawn,
%       with values corresponding to the local density of points, as
%       processed by the 'displayfunction'.
%       (default: false)
%     'displayfunction' <f/s>
%       When 'density' is true, a function that processes the density values
%       before displaying them. This can be any function handle, or one of
%       the following strings:
%         'atan':   arctangent, @atan
%         'linear': linear, same as the default, @(x) x
%         'log':    logarithmic, @(x) log(1 + x)
%         'sqrt':   square root, @(x) sqrt(1 + x)
%       Set to an empty matrix to mean no processing. Note that the density
%       values are counts (number of data points per pixel of output). Also
%       note that the function will be applied to the entire matrix of
%       densities.
%       (default: [])
%     'edge' <x>
%       The edge to leave around the data when 'tight' is true. This is
%       given in units of the range of the data in each axis.
%       (default: 0.1)
%     'labels' <b/v/c>
%       Labels to attach to the points, either as a numeric vector or as a
%       cell array of strings. If 'refseq' is not provided, the indices in
%       'labels' match those in x and y. If 'refseq' is provided, 'labels'
%       may just be set to true, to indicate that the labels from the
%       reference sequence should be used. Otherwise, if the reference
%       sequence mapping is towards integers, those integers will be used
%       as indices in the 'labels' array.
%       (default: no labels)
%     'maxpoints' <n>
%       Maximum number of points to display. If the number of elements in
%       the input vectors is larger than this, a uniform sampling will be
%       done so that at most n points will be plotted.
%       (default: 5000)
%     'range' <v>
%       Set the range of the figure, in the form [xmin, xmax, ymin, ymax].
%       (default: automatic)
%     'refseq' {r1, r2}
%       Assume that each of the vectors' components is mapped to a position
%       in a reference sequence by the arrays r1 and r2. In this case, the
%       lengths of x and y no longer need to be equal; instead, only those
%       components whose mapping to the reference sequence are common to
%       both vectors will be used for the scatterplot. In this case only,
%       the shape of the input matters: if the input is a matrix, then the
%       reference masks are applied to each dimension, instead of the
%       flattened data.
%       (default: no reference sequence mapping)
%     'resolution' <n>
%       Number of pixels to use for density plot (only when 'density' is
%       true). This can be a single number to use the same number of pixels
%       in each dimension, or a pair of numbers.
%       (default: 100)
%     'shapes' <c/v>
%       A single character representing the shape of the points in Matlab's
%       character format, or an array of characters giving the shape for
%       each point. This does not work when 'refseq' is provided.
%       (default: '.')
%     'sizes' <x/v>
%       A single number representing the size of the points, or a vector
%       giving the size of each of the points individually. This does not
%       work when 'refseq' is provided.
%       (default: 20)
%     'tight' <b>
%       Whether to use a framing that is tightly centered on the data, or
%       simply use Matlab's default.
%       (default: true)
%
% See also: PLOT.

% Tiberiu Tesileanu (2012-2015)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('colors', 'r', @(x) (ischar(x) && isvector(x)) || (isnumeric(x) && ismatrix(x)));
parser.addParameter('edge', 0.1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('labels', [], @(x) (isscalar(x) && islogical(x)) || (isvector(x) && (isnumeric(x) || iscell(x))));
parser.addParameter('maxpoints', 5000, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('range', [], @(x) isnumeric(x) && isvector(x) && length(x) == 4);
parser.addParameter('refseq', {}, @(x) iscell(x) && isvector(x) && length(x) == 2);
parser.addParameter('shapes', '.', @(x) isvector(x) && (ischar(x) || isnumeric(x)));
parser.addParameter('sizes', 20, @(x) isvector(x) && isnumeric(x));
parser.addParameter('tight', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('alpha', 1, @(x) isscalar(x) && isnumeric(x));

parser.addParameter('density', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('displayfunction', [], @(x) (ischar(x) && isvector(x)) || isa(x, 'function_handle'));
parser.addParameter('resolution', 100, @(n) isnumeric(n) && isscalar(n));

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle mapping to reference sequence
if ~isempty(params.refseq)
    [idxs1, idxs2] = findcommonidxs(params.refseq{1}, params.refseq{2});

    % treat matrices and vectors differently
    if ndims(x) == ndims(y) && ismatrix(x) && ~isvector(x)
        x = x(idxs1, idxs1);
        y = y(idxs2, idxs2);
    else
        x = x(idxs1);
        y = y(idxs2);
    end
    
    if ~isempty(params.labels)
        if islogical(params.labels)
            if params.labels
                params.labels = params.refseq{1}(idxs1);
            else
                params.labels = [];
            end
        else
            if isnumeric(params.refseq{1})
                params.labels = params.labels(params.refseq{1}(idxs1));
            else
                error([mfilename ':badlabels'], 'Can only mix ''labels'' with ''refseq'' when ''refseq'' is numeric.');
            end
        end
    end
else
    if islogical(params.labels)
        error([mfilename ':binlabel'], 'True or false ''labels'' is only allowed when ''refseq'' is provided.');
    end
end
if ~isempty(params.range)
    params.tight = false;
end
if isempty(params.displayfunction)
    params.displayfunction = @(x) x;
elseif ischar(params.displayfunction)
    switch params.displayfunction
        case 'linear'
            params.displayfunction = @(x) x;
        case 'log'
            params.displayfunction = @(x) log(1 + x);
        case 'sqrt'
            params.displayfunction = @(x) sqrt(1 + x);
        case 'atan'
            params.displayfunction = @atan;
        otherwise
            error([mfilename ':baddispf'], ['Unrecognized display function ' params.displayfunction '.']);
    end
end

% flatten the data
x = x(:);
y = y(:);

% some checks on the inputs
if length(x) ~= length(y)
    error([mfilename ':badlengths'], 'The sizes of x and y do not match.');
end
if ~isempty(params.labels) && length(params.labels) ~= length(x)
    error([mfilename ':badlabelno'], 'Number of labels doesn''t match number of points.');
end

if ischar(params.colors)
    ncolors = length(params.colors);
else
    ncolors = size(params.colors, 1);
end
nshapes = length(params.shapes);
nsizes = length(params.sizes);
if ncolors ~= 1 && ncolors ~= length(x)
    error([mfilename ':badcolors'], 'Bad size of ''colors'' option.');
end
if nshapes ~= 1 && nshapes ~= length(x)
    error([mfilename ':badshapes'], 'Bad size of ''shapes'' option.');
end
if nsizes ~= 1 && nsizes ~= length(x)
    error([mfilename ':badsizes'], 'Bad size of ''sizes'' option.');
end

if isscalar(params.resolution)
    params.resolution = [params.resolution params.resolution];
else
    if ~isvector(params.resolution) || length(params.resolution) ~= 2
        error([mfilename ':badres'], 'Resolution should be a vector of one or two nubmers.');
    end
end

% get the range of the data before subsampling
datarange = [min(x) max(x) min(y) max(y)];

% make sure we don't draw too many points
if ~params.density && length(x) > params.maxpoints
    step = ceil(length(x) / params.maxpoints);
    mask = (1:step:length(x));
    x = x(mask);
    y = y(mask);
    if ncolors > 1
        if ischar(params.colors)
            params.colors = params.colors(mask);
        else
            params.colors = params.colors(mask, :);
        end
        ncolors = length(mask);
    end
    if nshapes > 1
        params.shapes = params.shapes(mask);
        nshapes = length(mask);
    end
    if nsizes > 1
        params.sizes = params.sizes(mask);
        nsizes = length(mask);
    end
end

% figure out what the final range of the plot will be
if params.range
    final_range = params.range;
elseif params.tight
    edgex = params.edge*(datarange(2) - datarange(1));
    edgey = params.edge*(datarange(4) - datarange(3));
    final_range = datarange + [-edgex, edgex, -edgey, edgey];
end

xrange = (final_range(2) - final_range(1));
if xrange < eps
    xrange = 1;
    tmp = sum(final_range(1:2))/2;
    final_range(1:2) = tmp + [-xrange/2 xrange/2];
end
yrange = (final_range(4) - final_range(3));
if yrange < eps
    yrange = 1;
    tmp = sum(final_range(3:4))/2;
    final_range(3:4) = tmp + [-yrange/2 yrange/2];
end

if ~params.density
    % make the scatter plot
    if ncolors == 1 && nshapes == 1 && nsizes == 1 && params.alpha == 1
        plothandles = plot(x, y, ...
            'marker', params.shapes, ...
            'color', params.colors, ...
            'markersize', params.sizes, ...
            'linestyle', 'none');
    else
        % each plot is drawn in a different way
        if nshapes == 1
            params.shapes = repmat(params.shapes, 1, length(x));
        end
        if ncolors == 1
            if ischar(params.colors)
                params.colors = repmat(params.colors, 1, length(x));
            else
                params.colors = repmat(params.colors, length(x), 1);
            end
        end
        if nsizes == 1
            params.sizes = repmat(params.sizes, 1, length(x));
        end
        
        washold = ishold;
        hold on;
        
        plothandles = zeros(length(x), 1);
        for i = 1:length(x)
            if ischar(params.colors)
                crtcolor = params.colors(i);
            else
                crtcolor = params.colors(i, :);
            end
            if params.alpha == 1
                plothandles(i) = plot(x(i), y(i), ...
                    'marker', params.shapes(i), ...
                    'color', crtcolor, ...
                    'markersize', params.sizes(i), ...
                    'linestyle', 'none');
            else
                t = 0:pi/5:2*pi;
                shapex = params.sizes(i)*xrange/3000*cos(t);
                shapey = params.sizes(i)*yrange/3000*4/3*sin(t);
                plothandles(i) = patch(shapex + x(i), shapey + y(i), crtcolor, ...
                    'edgecolor', 'none');
                alpha(plothandles(i), params.alpha);
            end
        end
        
        if ~washold
            hold off;
        end
    end
else
    % make a density plot    
    xidx = min(1 + floor((x - final_range(1))*(params.resolution(1)/xrange)), params.resolution(1));
    yidx = min(1 + floor((y - final_range(3))*(params.resolution(2)/yrange)), params.resolution(2));
    densities = accumarray([yidx xidx], ones(size(xidx)), params.resolution);
    densities = params.displayfunction(densities);
    plothandles = image(final_range(1:2), final_range(3:4), densities, 'cdatamapping', 'scaled');
    set(gca, 'ydir', 'normal');
end

handles.plot = plothandles;
axis(final_range);

if ~isempty(params.labels)
    tmp = axis;
    range1 = tmp(2) - tmp(1);
    range2 = tmp(4) - tmp(3);
    shift1 = range1*0.01;
    shift2 = range2*0.01;
    for i = 1:length(params.labels)
        text(x(i) + shift1, y(i) + shift2, params.labels{i});
    end
end

end