function drawTernaryMixedBackground(varargin)
% drawTernaryMixedBackground Draw guiding axes through origin and guiding
% circles for pairs of ternary texture groups.
%   drawTernaryMixedBackground draws guiding axes going through the origin
%   and guiding circles for a plot representing a mix of two ternary
%   texture groups.
%
%   Options:
%    'circles'
%       Vector of radii at which to draw guiding circles.
%    'axes'
%       Set to false to not draw the coordinate axes.
%    'circleopts'
%       Plot options for the guiding circles.
%    'axisopts'
%       Plot options for the coordinate axes.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% parser.addParameter('circles', [1/4, 1/2, 3/4, 1], @(v) isempty(v) || (isvector(v) && isnumeric(v)));
parser.addParameter('circles', 1, @(v) isempty(v) || (isvector(v) && isnumeric(v)));
parser.addParameter('axes', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('circleopts', {'color', [0.8, 0.8, 0.8], 'linewidth', 0.5}, @(c) iscell(c));
parser.addParameter('axisopts', {'color', [0.8 0.8 0.8], 'linewidth', 0.5}, @(c) iscell(c));

% show defaults if asked
if nargin == 1 && strcmp(varargin{1}, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% we will overlay several graphical objects
wasHold = ishold;
hold on;

% we use this for drawing the circles
angleRange = linspace(0, 2*pi, 100);

% draw circles for orientation
circleRadii = params.circles;
for i = 1:length(circleRadii)
    radius = circleRadii(i);
    plot(radius*cos(angleRange), radius*sin(angleRange), params.circleopts{:});
end

% draw the main axes
plot([-1 1], [0 0], params.axisopts{:});
plot([0 0], [-1 1], params.axisopts{:});

% revert hold state
if ~wasHold
    hold off;
end

end
