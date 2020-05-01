function drawTernaryTriangle(varargin)
% drawTernaryTriangle Draw guiding simplex and circles for one of the
% ternary texture planes.
%   Options:
%    'edgelabels'
%       How to label the corners of the simplex. Set to
%           'probability'   to show probability vectors at simplex edges.
%           'digit'         to show ternary digit at simplex edges.
%           'none'          to skip labeling simplex edges.
%    'fontscale'
%       Scale factor for label fonts.
%    'innercircle'
%       Set to false to not draw the inner guiding circle.
%    'outercircle'
%       Set to false to not draw the outer guiding circle.
%    'simplex'
%       Set to false to not draw the simplex boundary.
%    'axes'
%       Set to false to not draw the coordinate axes.
%    'circleopts'
%       Plot options for guiding circles.
%    'simplexopts'
%       Plot options for simplex boundary.
%    'axisopts'
%       Plot options for coordinate axes.
%    'axisovershoot'
%       Fraction by which to overshoot the axes.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

[~, colorDict] = get_palette();

defaultColors = containers.Map;
defaultColors('circle') = lighten(colorDict('gray'), 0.75);
defaultColors('axis') = lighten(colorDict('gray'), 0.75);
defaultColors('simplex') = lighten(colorDict('dark blue'), 0.50);

parser.addParameter('edgelabels', 'probability', @(s) ismember(s, {'probability', 'digit', 'none'}));
parser.addParameter('fontscale', 1, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('innercircle', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('outercircle', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('simplex', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('axes', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('circleopts', {'color', defaultColors('circle'), 'linewidth', 0.5}, @(c) iscell(c));
parser.addParameter('simplexopts', {'color', defaultColors('simplex'), 'linewidth', 0.5}, @(c) iscell(c));
parser.addParameter('axisopts', {'color', defaultColors('axis'), 'linewidth', 0.5}, @(c) iscell(c));
parser.addParameter('axisovershoot', 0.5, @(x) isscalar(x) && isnumeric(x));

% show defaults if asked
if nargin == 1 && strcmp(varargin{1}, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% prepare to revert to the original hold state
washold = ishold;
hold on;

% half the size of a simplex edge
max_t2 = sqrt(3)/2;
% these will be used to draw the circles
angle_range = linspace(0, 2*pi, 100);

% draw circles for orientation (radius 1 and 1/2)
if params.outercircle
    plot(cos(angle_range), sin(angle_range), params.circleopts{:});
end
if params.innercircle
    plot(0.5*cos(angle_range), 0.5*sin(angle_range), params.circleopts{:});
end

% draw the probability triangle
if params.simplex
    plot([-1/2 1 -1/2 -1/2], [-max_t2 0 max_t2 -max_t2], params.simplexopts{:});
end

% draw the main axes
if params.axes
    f = 1 + params.axisovershoot;
    plot([0 f], [0 0], params.axisopts{:});
    plot([0 -1/2*f], [0 f*max_t2], params.axisopts{:});
    plot([0 -1/2*f], [0 -f*max_t2], params.axisopts{:});
end

% label the corners
if ~strcmp(params.edgelabels, 'none')
    % decide what kind of labels to use
    switch params.edgelabels
        case 'probability'
            l1 = '[0,1,0]';
            l2 = '[0,0,1]';
            l3 = '[1,0,0]';
        case 'digit'
            l1 = '1';
            l2 = '2';
            l3 = '0';
        otherwise
            error([mfilename ':badlabels'], 'Unknown edgelabels.');
    end
    
    % draw the labels
    h1 = text(1.05, 0, l1, 'fontsize', params.fontscale*12, 'color', [0.7 0.7 0.7]);
    h2 = text(-0.5,  max_t2, l2, 'fontsize', params.fontscale*12, 'color', [0.7 0.7 0.7]);
    h3 = text(-0.5, -max_t2, l3, 'fontsize', params.fontscale*12, 'color', [0.7 0.7 0.7]);
    
    % align the labels appropriately
    set(h1, 'horizontalalignment', 'left', 'verticalalignment', 'top');
    set(h2, 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
    set(h3, 'horizontalalignment', 'center', 'verticalalignment', 'top');
end

% restore hold state
if ~washold
    hold off;
end

end
