function [h_img, h_border] = bimage(varargin)
% BIMAGE Draw image with border.
%   BIMAGE(C) draws the array `C` as an image, placing the center of the
%   (i, j) pixel at position (j, i), and draws a border around it. See
%   IMAGE. The border is drawn using PLOT.
%
%   BIMAGE(x, y, C) specifies the image location, either in terms of the
%   (1, 1) corner, or the whole extent (if `x` and `y` are pairs). See
%   IMAGE.
%
%   [h_img, h_border] = BIMAGE(...) returns the handles to the image and
%   border that were created.
%
%   Name-value pairs are directly passed to IMAGE, except for the following:
%    'borderwidth'
%       Width of the border, passed to PLOT's 'linewidth'.
%       Set to 0 for no border.
%    'bordercolor'
%       Border color, as a string or 3-component RGB.
%    'gridwidth'
%       Width of grid lines between all elements; set to 0 for no lines.
%       This is directly passed to PLOT.
%    'gridcolor'
%       Grid color, as a string or 3-component RGB.
%
%   See also: IMAGE, PLOT.

% find first name-value pair, if it exists
idx = find(cellfun(@ischar, varargin), 1);
if ~isempty(idx)
    positional = varargin(1:idx-1);
    namevalues = varargin(idx:end);
else
    positional = varargin;
    namevalues = {};
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('borderwidth', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
% parser.addParameter('borderoverlap', 0.5, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('bordercolor', 'k', @(s) (ischar(s) && isscalar(s)) || isnumeric(x));
% parser.addParameter('borderalpha', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('gridwidth', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('gridcolor', 'k', @(s) (ischar(s) && isscalar(s)) || isnumeric(x));

if nargin == 1 && strcmp(varargin{1}, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(namevalues{:});
params = parser.Results;

% draw the image
h_img = image(positional{:}, parser.Unmatched);

% nothing more to do if we want no border and no grid
if params.borderwidth == 0 && params.gridwidth == 0
    h_border = [];
    return;
end

% figure out which axis we're using and turn hold on
ax = h_img.Parent;
was_hold = ishold(ax);
hold(ax, 'on');

% figure out the image extents on the axis
x_range = h_img.XData;
y_range = h_img.YData;
data = h_img.CData;
if isscalar(x_range)
    x_range = [x_range x_range + size(data, 2) - 1];
end
if isscalar(y_range)
    y_range = [y_range y_range + size(data, 1) - 1];
end

% make sure x_range and y_range are ordered
x_range = sort(x_range);
y_range = sort(y_range);

% correct the ranges to match pixel edges instead of centers
x_ps = (diff(x_range) + 1) / size(data, 2);
y_ps = (diff(y_range) + 1) / size(data, 1);
x_range = [x_range(1) - 0.5*x_ps x_range(2) + 0.5*x_ps];
y_range = [y_range(1) - 0.5*y_ps y_range(2) + 0.5*y_ps];

% draw the border
if params.borderwidth > 0
    h_border = plot(...
        [x_range(1) x_range(1) x_range(2) x_range(2) x_range(1)], ...
        [y_range(1) y_range(2) y_range(2) y_range(1) y_range(1)], ...
        'color', params.bordercolor, 'linewidth', params.borderwidth);
end

% draw the grid
if params.gridwidth > 0
    % draw horizontal lines
    yp = linspace(y_range(1), y_range(2), size(data, 1)+1);
    plot(repmat(x_range, size(data, 1)-1, 1)', repmat(yp(2:end-1)', 1, 2)', ...
        'color', params.gridcolor, 'linewidth', params.gridwidth);
    
    % draw vertical lines
    xp = linspace(x_range(1), x_range(2), size(data, 2)+1);
    plot(repmat(xp(2:end-1)', 1, 2)', repmat(y_range, size(data, 2)-1, 1)', ...
        'color', params.gridcolor, 'linewidth', params.gridwidth);
end

% revert hold state
if ~was_hold
    hold(ax, 'off');
end

end
