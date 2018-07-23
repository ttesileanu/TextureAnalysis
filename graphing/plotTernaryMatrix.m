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
%    'groupMaskFct'
%       A callable that takes in a group name and return true to keep the
%       group and false to dismiss it (and thus not draw it).
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
%       RGB colors.
%       Set to empty to disable per-subject coloring.
%    'limits'
%       The largest coordinate value to show in absolute value and in
%       either direction. This can be a pair of numbers, in which case the
%       first applies to single groups, and the second to mixed groups.
%    'plotterOptions'
%       Options to pass to `MatrixPlotter`.
%    'beautifyOptions'
%       Options to pass to beautifygraph.
%    Other options are sent directly to `plotTernaryThresholds`.
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
if nargin == 1 && strcmp(measurements, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% handle single set input
if ~iscell(measurements)
    measurements = {measurements};
end

% choose groups to display
if ~iscell(params.groupSelection)
    % collect group names according to the choice of groupSelection
    allGroups = {};
    for i = 1:length(measurements)
        % very basic input cleaning: remove whitespace from group names
        measurements{i}.groups = cellfun(@(s) s(~isspace(s)), ...
            measurements{i}.groups, 'uniform', false);
        
        % use the groupMaskFct to select the entries we want
        crtInitialGroups = unique(measurements{i}.groups);
        crtGroupMask = cellfun(params.groupMaskFct, crtInitialGroups);
        crtFinalGroups = crtInitialGroups(crtGroupMask);
        measurements{i} = selectMeasurements(measurements{i}, ...
            ismember(measurements{i}.groups, crtFinalGroups));
        
        % combine with allGroups
        if i == 1
            allGroups = crtFinalGroups;
        else
            switch params.groupSelection
                case 'intersect'
                    allGroups = intersect(allGroups, crtFinalGroups);
                case 'union'
                    allGroups = union(allGroups, crtFinalGroups);
                otherwise
                    error([mfilename ':badgrpsel'], 'Unknown option for groupSelection.');
            end
        end
    end
    uniqueGroups = sortGroups(unique(allGroups(:)));
else
    uniqueGroups = params.groupSelection;
end

% remove measurements entries that are not in uniqueGroups
for i = 1:length(measurements)
    measurements{i} = selectMeasurements(measurements{i}, ...
        ismember(measurements{i}.groups, uniqueGroups));
end

% figure out per-subject colors (if necessary)
if ~isempty(params.colorFct)
    if isa(params.colorFct, 'function_handle')
        % assign colors now, so that colors are consistent across different
        % groups even if some subjects don't have data for some groups
        allSubjects = {};
        for i = 1:length(measurements)
            if isfield(measurements{i}, 'subjects')
                allSubjects = union(allSubjects, measurements{i}.subjects);
            end
        end
        allSubjectColors = params.colorFct(length(allSubjects));
        colorFct = containers.Map(allSubjects, num2cell(allSubjectColors, 2));
    else
        colorFct = params.colorFct;
    end
else
    colorFct = [];
end

% use the MatrixPlotter to make a figure containing a matrix of plots
plotter = MatrixPlotter(length(uniqueGroups), params.plotterOptions{:});
plotter.axAspect = 1;
while plotter.next
    % find the index of the current plot
    i = plotter.index;
    
    % plot the measurements for this plane
    hold on;
    for k = 1:length(measurements)
        crtMeasurements = selectMeasurements(measurements{k}, ...
            strcmp(measurements{k}.groups, uniqueGroups{i}));
        if isfield(crtMeasurements, 'thresholdIntervals')
            errors = crtMeasurements.thresholdIntervals;
        else
            errors = [];
        end
        plotTernaryThresholds(crtMeasurements, ...
            'marker', selectMod(params.markers, k), ...
            'color', selectMod(params.colors, k), ...
            'size', selectMod(params.sizes, k), ...
            'errors', errors, ...
            'colorFct', colorFct, ...
            'drawGuides', k == 1, ...
            unmatched);
    end
    
    % set axis limits and beautify
    nGroups = 1 + sum(uniqueGroups{i} == ';');
    switch nGroups
        case 1
            xlim([-params.limits(1) params.limits(1)]);
            ylim([-params.limits(1) params.limits(1)]);
        case 2
            xlim([-params.limits(2) params.limits(2)]);
            ylim([-params.limits(2) params.limits(2)]);
    end

    beautifygraph(params.beautifyOptions{:});
end

% set up figure
fig = plotter.getFigure;
fig.Name = [fig.Name ' (crosses = PP, dots = NI)'];

preparegraph;

end