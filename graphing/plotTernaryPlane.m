function plotTernaryPlane(meas, mask, varargin)
% plotTernaryPlane Plot thresholds in one ternary texture plane (simple or 
% mixed).
%   plotTernaryPlane(meas, mask) plots the measurements from the structure
%   `meas` corresponding to the `true` entries in the logical array `mask`
%   on an appropriate background depending on whether the measurements are
%   for single or mixed planes. The `meas` structure needs at least fields
%   `directions` and `thresholds`. Error bars are drawn if the field
%   `thresholdIntervals` also exists. Per-subject coloring can be used if
%   the field `subjects` is present. For mixed-planes, a `groups` field
%   must also be present in `meas` to properly label the axes. For both
%   single and mixed planes, a `groups` field must be present to set a
%   title.
%
%   Options:
%    'symbol'
%       Marker (symbol) to use.
%    'symbolColor'
%       3-component RGB vector giving the color of the symbols to use.
%    'symbolSize'
%       Symbol size to use.
%    'errorBars'
%       Set to true to draw error bars, false otherwise.
%    'errorColor'
%       3-component RGB vector giving the color of the error bars.
%    'ellipses'
%       Set to true to show best-fit ellipses, false otherwise.
%    'colorFct'
%       Colormap function giving the colors to use when there are several
%       subjects. This can be one of the following:
%         (a) a colormap-like function, such that `colorFct(n)` returns an
%             n x 3 matrix of RGB colors.
%         (b) an instance of `containers.Map`, mapping subject names to
%             RGB colors.
%       Set to empty to disable per-subject coloring.
%    'beautifyOptions'
%       Options to pass to beautifygraph.
%    'triangleOptions'
%       Options to pass to drawTernaryTriangle.
%    'limits'
%       The largest coordinate value to show in absolute value and in
%       either direction. This can be a pair of numbers, in which case the
%       first applies to single groups, and the second to mixed groups.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('symbol', 'x', @(s) ischar(s));
parser.addParameter('symbolColor', [0.8 0.3 0.3], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('symbolSize', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('errorBars', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('errorColor', [0.5 0.5 0.5], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('ellipses', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('colorFct', @parula, @(f) isempty(f) || isa(f, 'function_handle') || isa(f, 'containers.Map'));
parser.addParameter('beautifyOptions', {'box', 'on', 'tickdir', 'in', 'titlesize', 12}, ...
    @(c) iscell(c) && (isempty(c) || isvector(c)));
parser.addParameter('limits', [2 1], @(v) isnumeric(v) && isvector(v) && ismember(length(v), [1 2]));
parser.addParameter('triangleOptions', {}, @(c) iscell(c) && isvector(c));

% show defaults if requested
if nargin == 1 && strcmp(meas, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% handle 1-element limits
if length(params.limits) == 1
    params.limits = repmat(params.limits, 1, 2);
end

% select the elements that match the mask
maskedDirections = meas.directions(mask);
maskedThresholds = meas.thresholds(mask);
if isfield(meas, 'thresholdIntervals') && params.errorBars
    maskedThresholdIntervals = meas.thresholdIntervals(mask, :);
else
    maskedThresholdIntervals = [];
end
if isfield(meas, 'subjects')
    maskedSubjects = meas.subjects(mask);
else
    maskedSubjects = [];
end

% figure out how whether this is a single-group or a mixed-group plane
if mod(maskedDirections{1}, 3) ~= 0
    error([mfilename ':baddirsz'], 'Directions vector sizes aren''t compatible with ternary planes.');
end
nGroups = length(maskedDirections{1}) / 3;
if any(cellfun(@(c) length(c) ~= nGroups*3, maskedDirections))
    error([mfilename ':mixeddirs'], 'Not all directions have the same size.');
end

% check that this is a supported number of groups, and draw the appropriate
% background
switch nGroups
    case 1
        drawTernaryTriangle(params.triangleOptions{:});
    case 2
        if ~isfield(meas, 'groups')
            groupOpts = {'Group 1', 'Group 2'};
        else
            maskedGroups = meas.groups(mask);
            groupOpts = strsplit(maskedGroups{1}, ';');
        end
        drawTernaryMixedBackground(groupOpts{:});
    otherwise
        error([mfilename ':badngrp'], 'This function only works with single groups and pairs of groups.');
end

% preserve hold state, to reset later
wasHold = ishold;
hold on;

% make a dictionary of subject colors, if we have subject info
havePerSubjectColor = false;
if ~isempty(maskedSubjects)
    allSubjects = unique(maskedSubjects);
else
    allSubjects = 'subject';
end
if length(allSubjects) > 1
    if isempty(params.colorFct)
        % same color for all subjects
        colorMap = containers.Map(allSubjects, ...
            repmat({params.symbolColor}, 1, length(allSubjects)));
    elseif isa(params.colorFct, 'function_handle')
        % pick out colors from a color matrix
        allSubjectColors = params.colorFct(length(allSubjects));
        colorMap = containers.Map(allSubjects, num2cell(allSubjectColors, 2));
        havePerSubjectColor = true;
    else
        % the user already provided a dictionary
        colorMap = params.colorFct;
        havePerSubjectColor = true;
    end
else
    colorMap = containers.Map(allSubjects, {params.symbolColor});
end

% draw the thresholds
for i = 1:length(allSubjects)
    % select the measurements for one subject
    subject = allSubjects{i};
    if ~isempty(maskedSubjects)
        subjectMask = strcmp(maskedSubjects, subject);
        subjectThresholds = maskedThresholds(subjectMask);
        subjectDirections = maskedDirections(subjectMask);
        if ~isempty(maskedThresholdIntervals)
            subjectThresholdIntervals = maskedThresholdIntervals(subjectMask, :);
        else
            subjectThresholdIntervals = [];
        end
    else
        subjectThresholds = maskedThresholds;
        subjectDirections = maskedDirections;
        subjectThresholdIntervals = maskedThresholdIntervals;
    end
    
    % draw
    subjectColor = colorMap(subject);
    if havePerSubjectColor
        errorColor = mixcolor(subjectColor, params.errorColor);
    else
        errorColor = params.errorColor;
    end
    plotTernaryThresholds(subjectThresholds, subjectDirections, ...
        'errors', subjectThresholdIntervals, 'color', subjectColor, ...
        'size', params.symbolSize, 'errorColor', errorColor, ...
        'ellipse', params.ellipses);
end

% arrange plot sizes and labels
axis equal;

switch nGroups
    case 1
        xlim([-params.limits(1) params.limits(1)]);
        ylim([-params.limits(1) params.limits(1)]);
    case 2
        xlim([-params.limits(2) params.limits(2)]);
        ylim([-params.limits(2) params.limits(2)]);
end

% reset hold state
if ~wasHold
	hold off;
end

% set title, beautify
if isfield(meas, 'groups')
    maskedGroups = meas.groups(mask);
    title(maskedGroups{1});
end

beautifygraph(params.beautifyOptions{:});

end