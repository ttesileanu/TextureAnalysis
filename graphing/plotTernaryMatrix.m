function [plotter, uniqueGroups] = plotTernaryMatrix(measurements, varargin)
% plotTernaryMatrix Make plots showing several sets of ternary texture
% measurements.
%   plotTernaryMatrix(measurements) makes a plot for each texture group
%   present in the `measurements`. This structure should contain fields
%   `groups`, `directions`, and `thresholds`. If `thresholdIntervals` is
%   also present, it is used to draw error bars (unless error bars are
%   disabled; see options below).
%
%   plotTernaryMatrix({measurements1, measurements2, ...}) draws several
%   sets of measurements. By default, all groups that appear in at least
%   one measurement set are displayed, but see 'groupSelection' option
%   below.
%
%   [plotter, uniqueGroups] = plotTernaryMatrix(...) returns the
%   `MatrixPlotter` instance used for plotting and the list of groups
%   mapping to the list of axes in the plotter.
%
%   By default the plots are arranged in a rectangular grid that covers
%   most of the screen area. This can be changed using the 'plotterOptions'
%   below.
%
%   Options:
%    'groupSelection'
%       Set either to a cell array of group names -- in which case these
%       are the only groups that will be drawn. Or set to
%         'intersect':  draw only groups that appear in all measurement sets
%         'union':      draw all groups that appear in at least one set
%    'colors'
%       Cell array of colors (as 3-component RGB vectors) to be used to
%       display the measurements. These are used cyclically if there are
%       more measurements than colors.
%    'markers'
%       Cell array of markers to be used for the measurements. These are
%       used cyclically if there are more measurements than markers.
%    'sizes'
%       Vector of marker sizes. These are used cyclically if there are more
%       measurements than sizes.
%    'colorFct'
%       Colormap function giving the colors to use when there are several
%       subjects. This overrides the 'colors' option. It can be one of the
%       following:
%         (a) a colormap-like function, such that `colorFct(n)` returns an
%             n x 3 matrix of RGB colors.
%         (b) an instance of `containers.Map`, mapping subject names to
%             RGB colors.
%       Set to empty to disable per-subject coloring.
%    Other options are sent to `showTernaryChange`.
%
%   See also: plotTernaryThresholds.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('groupMaskFct', @(g) true);
parser.addParameter('groupSelection', 'union', @(c) (iscell(c) && isvector(c)) ...
    || ismember(c, {'intersect', 'union'}));
parser.addParameter('colors', {[0 0.3438 0.7410], [0.8 0.3 0.3]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('markers', {'.', 'x'}, @(c) iscell(c) && isvector(c));
parser.addParameter('sizes', [8, 5], @(x) isnumeric(x) && isvector(x) && all(x > 0));
parser.addParameter('colorFct', [], @(f) isempty(f) || isa(f, 'function_handle') ...
    || isa(f, 'containers.Map'));
parser.addParameter('limits', [2 1], @(v) isvector(v) && isnumeric(v) && ...
    ismember(numel(v), [1 2]));
parser.addParameter('plotterOptions', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('beautifyOptions', {}, @(c) iscell(c) && isvector(c));

% show defaults if requested
if nargin == 1 && isequal(measurements, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% process options
[plotter, uniqueGroups] = showTernaryChange(...
    measurements, {}, 'groupMaskFct', params.groupMaskFct, ...
    'groupSelection', strcat(params.groupSelection, 'All'), ...
    'colorsBefore', params.colors, 'markersBefore', params.markers, ...
    'sizesBefore', params.sizes, 'colorFctBefore', params.colorFct, ...
    'limits', params.limits, 'plotterOptions', params.plotterOptions, ...
    'beautifyOptions', params.beautifyOptions, unmatched);

end
