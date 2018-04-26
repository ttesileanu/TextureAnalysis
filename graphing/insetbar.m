function insetbar(position, heights, varargin)
% INSETBAR Overlay bar in an existing plot.
%   INSETBAR(position, heights) draws a bar plot on the current axes,
%   starting at `position(1:2)` in data coordinates, and extending for
%   `position(3:4)`, again in data coordinates. The `heights` are
%   normalized so that a value of 1 touches the upper limit of the inset.
%
%   INSETBAR(position, heights, width) sets the width of the bars, as a
%   fraction of the spacing between bars.
%
%   XXX This currently gives undesirable results if some of the heights are
%   negative.
%
%   Options:
%    'color'
%       Color of the bars, as a character, or as RGB vector.
%    'bottomcolor'
%       Color of the edge drawn at the bottom of the plot, as a character
%       or RGB vector. If empty, then no edge is drawn.
%    'bottomdist'
%       Distance from bottom edge of the bar plot to the bottom of the
%       bars, in units of the total height of the bar plot.
%    'edge'
%       Set to true to draw edges on the bars, false otherwise.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

check_color = @(c) isempty(c) || (ischar(c) && length(c) == 1) || ...
    (isnumeric(c) && isvector(c) && length(c) == 3);
parser.addOptional('width', 0.8, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('color', [0.5 0.7 1], check_color);
parser.addParameter('bottomcolor', 'k', check_color);
parser.addParameter('bottomdist', 0.1, @(x) isscalar(x) && isnumeric(x) && x >= 0);
parser.addParameter('edge', false, @(b) islogical(b) && isscalar(b));

if nargin == 1 && strcmp(varargin{1}, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% draw
n = length(heights);
if n > 0
    patch_opts = {};
    if ~params.edge
        patch_opts = [patch_opts {'edgecolor', 'none'}];
    end
    
    w0 = position(3)/n;
    centers = w0/2 + (0:n-1)*w0;
    
    % XXX not the most efficient...
    w = params.width*w0;
    y0 = position(2) + params.bottomdist*position(4);
    h0 = (1 - params.bottomdist)*position(4);
    scaling = h0;
    for i = 1:length(heights)
        x0 = position(1) + centers(i) - w/2;
        h = heights(i)*scaling + eps;
        fill([x0 x0 x0+w x0+w x0], [y0 y0+h y0+h y0 y0], params.color, patch_opts{:});
    end
end

% draw edge at bottom, if necessary
if ~isempty(params.bottomcolor)
    plot([position(1) position(1) + position(3)], repmat(position(2), 1, 2), ...
        'color', params.bottomcolor);
end

end