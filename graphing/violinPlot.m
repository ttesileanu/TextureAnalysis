function h = violinPlot(x, y, varargin)
% violinPlot Make a violin plot.
%   violinPlot(x, y) plots a kernel-density estimate of the `y`
%   distribution at each value of the categorical variable `x`. When `x` is
%   not numeric, all values of `x` are sorted and assinged locations on the
%   x-axis of the plot from 1 to `length(unique(x))`. The different
%   categories are shown in different colors according to these cateogries.
%   The colors are by default drawn from the 'lines' colormap (but see
%   options).
%
%   h = violinPlot(x, y) returns the handles for the plotted data.
%
%   Options:
%     'ax'
%       Axes to draw in.
%     'colors'
%       Colormap to use, either as a string, or as a 3-column RGB matrix.
%       Colors are chosen sequentially from the colormap, cycling back to
%       the beginning if necessary.
%     'alpha'
%       Opacity of violins.
%     'width'
%       Width of each violin.
%     'boxes'
%       Set to `true` to draw interquartile ranges and medians.
%     'boxOpts'
%       Options to pass to plot for the interquartile range boxes and
%       medians.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('ax', []);
parser.addParameter('colors', 'lines', @(c) (ischar(c) && isvector(c)) || ...
    (ismatrix(c) && size(c, 2) == 3 && isnumeric(c) && isreal(c)));
parser.addParameter('alpha', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('width', 0.5, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('boxes', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('boxOpts', {'linewidth', 0.5, 'color', 'k'}, @(c) iscell(c));

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
        error([mfilename ':badcmap'], 'Invalid colormap name passed to colors option.');
    end
    
    params.colors = eval(params.colors);
end

% calculate the kernel density estimates
categories = unique(x);
kdes = cell(1, length(categories));
c = zeros(length(categories), 3);
category_xs = zeros(length(categories), 1);
for i = 1:length(categories)
    if ~isnumeric(x)
        category_xs(i) = i;
    else
        category_xs(i) = categories(i);
    end
    
    mask = ismember(x, categories(i));
    allY = y(mask);
    [kde_x, kde_y] = ksdensity(allY);
    kde_x = 0.5 * params.width * kde_x / max(kde_x);
    kdes{i} = {kde_x(:), kde_y(:)};
    c(i, :) = params.colors(mod(i - 1, size(params.colors, 1)) + 1, :);
end

% draw violins
h = [];
for i = 1:length(kdes)
    crt_kde = kdes{i};
    crt_h = fill(category_xs(i) + [crt_kde{1} ; -flipud(crt_kde{1})], ...
        [crt_kde{2} ; flipud(crt_kde{2})], ...
        c(i, :), 'facealpha', params.alpha, 'edgecolor', 0.5*c(i, :));
    h = [h crt_h]; %#ok<AGROW>
end

% draw boxes, if necessary
if params.boxes
    wasHold = ishold;
    hold on;
    
    % decide how wide the boxes should be
    boxWidth = params.width;
    
    for i = 1:length(categories)
        mask = ismember(x, categories(i));
        allY = y(mask);
        crt_x = category_xs(i);
        
        quartiles = quantile(allY, [0.25, 0.50, 0.75]);
        
        % draw the box
        plot(crt_x + 0.5*boxWidth * [-1 1 1 -1 -1], ...
             [quartiles(1) quartiles(1) quartiles(3) quartiles(3) quartiles(1)], ...
             params.boxOpts{:});
         
        % draw the median
        plot(crt_x + 0.5*boxWidth * [-1 1], ...
             quartiles(2)*[1 1], params.boxOpts{:});
    end
    
    if ~wasHold
        hold off;
    end
end

end
