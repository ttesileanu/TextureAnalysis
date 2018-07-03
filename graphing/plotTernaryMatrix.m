function plotTernaryMatrix(pred, meas, varargin)
% plotTernaryMatrix Make plots for each texture group comparing predicted
% to measured thresholds.
%   plotTernaryMatrix(pred, meas) makes a plot for each texture group
%   present in the measurement structure `meas`. This structure should
%   contain `groups`, `directions`, and `thresholds`. If `thresholdIntervals`
%   are present, they are also used (unless error bars are disabled; see
%   options below). The predicted thresholds `pred` should be given as a
%   vector with the same ordering as that in the measurements structure.
%
%   plotTernaryMatrix([], meas) plots only the measurements.
%
%   By default the plots are arranged in a rectangular grid that covers
%   most of the screen area. This can be changed using the 'plotterOpts'
%   below.
%
%   Options:
%    'mask'
%       Boolean mask showing which measurements to include in the plots.
%    'predColor'
%       3-component RGB vector giving the color of the predictions.
%    'measColor'
%       3-component RGB vector giving the color of the measurements when
%       'colorFct' isn't used.
%    'measColorFct'
%       Colormap function giving the colors to use when there are several
%       subjects. This can be one of the following:
%         (a) a colormap-like function, such that `colorFct(n)` returns an
%             n x 3 matrix of RGB colors.
%         (b) an instance of `containers.Map`, mapping subject names to
%             RGB colors.
%       Set to empty to disable per-subject coloring.
%    'predSize'
%       Symbol size for predictions.
%    'measSize'
%       Symbol size for measurements.
%    'multi'
%       How to deal with mixed groups. If set to `false`, mixed groups are
%       ignored. If set to `true`, all groups are drawn. It can also be a
%       vector giving the multiplicities to draw, for instance [1] shows
%       only single groups, [2] shows pairs, and [1, 2] shows both. Note
%       that (for now?) only single groups and pairs are supported.
%    'plotterOpts'
%       Options to pass to MatrixPlotter.
%   Other options are sent directly to plotTernaryPlane.
%
%   See also: plotTernaryPlane.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

% check for 'defaults' request
if nargin == 1 && strcmp(pred, 'defaults')
    showDefaults = true;
    % some of the defaults below depend on the amount of measurement data,
    % so we need to provide a mock `meas` structure
    meas.thresholds = [];
else
    showDefaults = false;
end

parser.addParameter('mask', true(size(meas.thresholds)), @(b) islogical(b) && isvector(b) && length(b) == length(meas.thresholds));
parser.addParameter('predColor', [0 0.3438 0.7410], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('measColor', [0.8 0.3 0.3], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('measColorFct', @parula, @(f) isempty(f) || isa(f, 'function_handle') || isa(f, 'containers.Map'));
parser.addParameter('predSize', 8, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('measSize', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('multi', false, @(b) (islogical(b) && isscalar(b)) | ...
    (isnumeric(b) && isvector(b)));
parser.addParameter('plotterOpts', {}, @(c) iscell(c) && isvector(c));

% show defaults if requested
if showDefaults
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% figure out what to do with multiple groups
if islogical(params.multi)
    if params.multi
        params.multi = [1, 2];
    else
        params.multi = 1;
    end
end

% figure out which entries are tuples, and of how many elements
multiplicities = 1 + cellfun(@(s) sum(s == ';'), meas.groups);
% figure out what to keep
multiMask = ismember(multiplicities, params.multi);

% very basic input cleaning: remove whitespace from group names
meas.groups = cellfun(@(s) s(~isspace(s)), meas.groups, 'uniform', false);

% figure out which groups we're showing
mask = params.mask & multiMask;
maskedGroups = meas.groups(mask);
% these are all the groups we have; sort them by order, then alphabetically
uniqueGroups = sortGroups(unique(maskedGroups(:)));

% figure out per-subject colors (if necessary)
if ~isempty(params.measColorFct) && isfield(meas, 'subjects')
    if isa(params.measColorFct, 'function_handle')
        % assign colors now, so that colors are consistent across different
        % groups even if some subjects don't have data for some groups
        % pick out colors from a color matrix
        allSubjects = unique(meas.subjects(params.mask & multiMask));
        allSubjectColors = params.measColorFct(length(allSubjects));
        colorFct = containers.Map(allSubjects, num2cell(allSubjectColors, 2));
    else
        colorFct = params.measColorFct;
    end
else
    colorFct = [];
end

% use the MatrixPlotter to make a figure containing a matrix of plots
plotter = MatrixPlotter(length(uniqueGroups), params.plotterOpts{:});
plotter.axAspect = 1;
while plotter.next
    % find the index of the current plot
    i = plotter.index;
    
    % plot the measurements for this plane
    displayMask = mask & strcmp(meas.groups, uniqueGroups{i});
    plotTernaryPlane(meas, displayMask, 'symbolColor', params.measColor, ...
        'symbolSize', params.measSize, unmatched);
    
    % plot the predictions, if any
    if ~isempty(pred)
        hold on;
        plotTernaryThresholds(pred(displayMask), meas.directions(displayMask), ...
            'color', params.predColor, 'colorFct', colorFct);
    end
end

fig = plotter.getFigure;
fig.Name = [fig.Name ' (crosses, open circles = PP, dots = NI)'];

preparegraph;

end