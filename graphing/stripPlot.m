function h = stripPlot(x, y, varargin)
% stripPlot Scatter plot where one variable is categorical.
%   stripPlot(x, y) makes a scatterplot of `y` vs. `x` in which `x` is
%   assumed categorical. When `x` is not numeric, all values of `x` are
%   sorted and assinged locations on the x-axis of the plot from 1 to
%   `length(unique(x))`. The different categories are shown in different
%   colors according to these cateogries. The colors are by default drawn
%   from the 'lines' colormap (but see options). Jitter can be added to
%   separate the points along the categorical direction (see options).
%
%   h = stripPlot(x, y) returns the handles for the plotted data.
%
%   Options:
%     'ax'
%       Axes to draw in.
%     'colors'
%       Colormap to use, either as a string, or as a 3-column RGB matrix.
%       Colors are chosen sequentially from the colormap, cycling back to
%       the beginning if necessary.
%     'sizes'
%       Point sizes. XXX For now this should be a scalar.
%     'jitter'
%       Amount of jitter to use.
%     'kde'
%       Scale amount of jitter by a kernel density estimate of the data in
%       each category.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('ax', []);
parser.addParameter('colors', 'lines', @(c) (ischar(c) && isvector(c)) || ...
    (ismatrix(c) && size(c, 2) == 3 && isnumeric(c) && isreal(c)));
parser.addParameter('sizes', [], @(s) isempty(s) || (isnumeric(s) && isscalar(s)));
parser.addParameter('jitter', 0, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('kde', false, @(b) islogical(b) && isscalar(b));

% show defaults if requested
if nargin == 1 && strcmp(x, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% handle string colormap
if ischar(params.colors)
    % get the colors from the colormap
    % !!! we need to use `eval`; for some safety, check that the name given
    % is a simple variable/function name
    if ~isvarname(params.colors)
        error([mfilename ':badcmap'], 'Invalid colormap name passeed to colors option.');
    end
    
    params.colors = eval(params.colors);
end

% create a color matrix
categories = unique(x);
c = zeros(length(x), 3);
xAll = zeros(size(y));
for i = 1:length(categories)
    mask = ismember(x, categories(i));
    if ~isnumeric(x)
        xAll(mask) = i;
    else
        xAll(mask) = categories(i);
    end
    color = params.colors(mod(i-1, size(params.colors, 1)) + 1, :);
    for k = 1:3
        % using `ismember` allows this to work seemlessly with either numeric
        % or string data
        c(mask, k) = color(k);
    end
end

% add jitter
if ~params.kde
    xAll = xAll + params.jitter*(rand(size(xAll)) - 0.5);
else
    for i = 1:length(categories)
        mask = ismember(x, categories(i));
        allY = y(mask);
        kde = ksdensity(allY, allY);
        kde = kde / max(kde);
        crtJitter = params.jitter*(rand(sum(mask), 1) - 0.5) .* kde;
        xAll(mask) = xAll(mask) + crtJitter;
    end
end

% draw
h = scatter(xAll, y, params.sizes, c);

end