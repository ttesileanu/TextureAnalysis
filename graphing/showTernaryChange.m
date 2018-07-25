function showTernaryChange(measurementsBefore, measurementsAfter, varargin)
% showTernaryChange Show a change in ternary threshold measurements.
%   showTernaryChange(measurementsBefore, measurementsAfter) makes a plot
%   of `measurementsBefore` and `measurementsAfter` (similar to
%   `plotTernaryMatrix`), and adds lines connecting the 'before' and
%   'after' measurements.
%
%   In general, the groups in `measurementsBefore` and `measurementsAfter`
%   may differ. Groups that are present only in `measurementsAfter` are not
%   shown (but see the 'groupSelection' option). Measurements that don't
%   have an 'after' are drawn in a different color (see options below).
%
%   Several sets of 'before' and 'after' measurements can be shown by
%   passing cell arrays to `measurementsBefore` and `measurementsAfter`. In
%   this case, the number of elements in the two arrays should be the same.
%
%   Options:
%    'groupMaskFct'
%       A callable that takes in a group name and returns true to keep the
%       group and false to dismiss it (and thus not draw it).
%    'groupSelection'
%       Set either to a cell array of group names -- in which case these
%       are the only groups that will be drawn. Or set to
%        'intersectBefore': 
%           show only groups that appear in all 'before' measurement sets
%        'unionBefore':
%           show all groups that appear in at least one 'before' set
%        'intersectAfter': 
%           show only groups that appear in all 'after' measurement sets
%        'unionAfter':
%           show all groups that appear in at least one 'after' set
%        'intersectAll': 
%           show only groups that appear in all measurement sets, 'before'
%           and 'after'
%        'unionAll':
%           show all groups that appear in at least one set, either
%           'before' or 'after'
%    'colorsBefore'
%    'colorsAfter'
%       Cell array of colors (as 3-component RGB vectors) to be used to
%       display the 'before' and 'after' measurements, respectively. These
%       are used cyclically if there are more measurements than colors.
%    'colorsInvariant'
%       Cell array of colors to use for each measurement set for those
%       thresholds that were invariant under the transformation that took
%       place between 'before' and 'after'. If the 'after' measurements
%       have a `changed` field, that field is used to determine which
%       thresholds were not invariant. Otherwise, measurements are assume
%       to be invariant if they are very close (see 'invarianceTol' below).
%    'colorsNoMatch'
%       Cell array of colors to use for each measurement set for those
%       thresholds that have no match between 'before' and 'after'.
%    'lineColors'
%       Cell array of colors to use to connect 'before' and 'after'
%       measurements.
%    'markersBefore'
%    'markersAfter'
%       Cell array of markers to be used for the 'before' and 'after'
%       measurements, respectively. These are used cyclically if there are
%       more measurements than markers.
%    'sizesBefore'
%    'sizesAfter'
%       Vector of marker sizes for 'before' and 'after' measurements',
%       respectively. These are used cyclically if there are more
%       measurements than sizes.
%    'colorFctBefore'
%    'colorFctAfter'
%       Colormap functions giving the colors to use when there are several
%       subjects, for 'before' and 'after' measurements, respectively.
%       These override the 'colorsBefore'/'colorsAfter' options, and they
%       should be callables that take in the name of a subject and return
%       an RGB color. Set to empty to disable per-subject coloring.
%    'invarianceTol'
%       How similar two thresholds have to be to be considered identical.
%       This is used to determined invariant directions when the 'after'
%       measurements don't contain a `changed` field.
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
%   See also: MatrixPlotter, plotTernaryMatrix.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('groupMaskFct', @(g) true);
parser.addParameter('groupSelection', 'unionBefore', ...
    @(c) (iscell(c) && isvector(c)) || ismember(c, ...
    {'intersectBefore', 'unionBefore', 'intersectAfter', 'unionAfter', ...
     'intersectAll', 'unionAll'}));
parser.addParameter('colorsBefore', {[0.3 0.6438 1.0000], [1.0 0.6 0.6]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('colorsAfter', {[0 0.3438 0.7410], [0.8 0.3 0.3]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('colorsInvariant', {[0.9 0.9 0.9]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('colorsNoMatch', {[0.2 0.2 0.2]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('lineColors', {[0 0.4938 0.8910], [0.95 0.45 0.45]}, ...
    @(c) iscell(c) && isvector(c));
parser.addParameter('markersBefore', {'.', 'x'}, @(c) iscell(c) && isvector(c));
parser.addParameter('markersAfter', {'.', 'x'}, @(c) iscell(c) && isvector(c));
parser.addParameter('sizesBefore', [8, 5], @(x) isnumeric(x) && isvector(x) && all(x > 0));
parser.addParameter('sizesAfter', [8, 5], @(x) isnumeric(x) && isvector(x) && all(x > 0));
parser.addParameter('colorFctBefore', []);
parser.addParameter('colorFctAfter', []);
parser.addParameter('invarianceTol', 1e-8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('limits', [2 1], @(v) isvector(v) && isnumeric(v) && ...
    ismember(numel(v), [1 2]));
parser.addParameter('plotterOptions', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('beautifyOptions', {}, @(c) iscell(c) && isvector(c));

% show defaults if requested
if nargin == 1 && strcmp(measurementsBefore, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% XXX there is a lot of code repetition here from plotTernaryMatrix
% XXX probably the best way to fix this is to implement plotTernaryMatrix
%     as a special case of this function, with no 'after' measurements

% handle single set inputs
if ~iscell(measurementsBefore)
    measurementsBefore = {measurementsBefore};
end
if ~iscell(measurementsAfter)
    measurementsBefore = {measurementsAfter};
end

if length(measurementsBefore) ~= length(measurementsAfter)
    error([mfilename ':lenmismatch'], 'Should have matching numbers of  ''before'' and ''after'' measurements.');
end

% very basic input cleaning: remove whitespace from group names
for i = 1:length(measurementsBefore)
    measurementsBefore{i}.groups = cellfun(@(s) s(~isspace(s)), ...
        measurementsBefore{i}.groups, 'uniform', false);
end
for i = 1:length(measurementsAfter)
    measurementsAfter{i}.groups = cellfun(@(s) s(~isspace(s)), ...
        measurementsAfter{i}.groups, 'uniform', false);
end

% choose groups to display
if ~iscell(params.groupSelection)
    % collect group names according to the choice of groupSelection
    % first figure out whether we need to look only at 'before', only at
    % 'after', or at everything, and separate that from 'intersect' vs.
    % 'union'
    switch params.groupSelection
        case {'intersectBefore', 'unionBefore'}
            measForSel = measurementsBefore;
            groupSelShort = params.groupSelection(1:end-6);
        case {'intersectAfter', 'unionAfter'}
            measForSel = measurementsAfter;
            groupSelShort = params.groupSelection(1:end-5);
        case {'intersectAll', 'unionAll'}
            measForSel = [measurementsBefore measurementsAfter];
            groupSelShort = params.groupSelection(1:end-3);
        otherwise
            error([mfilename ':badgrpsel'], 'Unknown option for groupSelection.');
    end
    % now collect group names
    allGroups = {};
    for i = 1:length(measForSel)
        crtGroups = measForSel{i}.groups;
        
        % combine with allGroups
        if i == 1
            allGroups = crtGroups;
        else
            switch groupSelShort
                case 'intersect'
                    allGroups = intersect(allGroups, crtGroups);
                case 'union'
                    allGroups = union(allGroups, crtGroups);
                otherwise
                    error([mfilename ':badgrpsel'], 'Unknown option for groupSelection.');
            end
        end
    end
    
    % select groups according to the groupMaskFct
    allGroupsMask = cellfun(params.groupMaskFct, allGroups);
    allGroupsFinal = allGroups(allGroupsMask);
    
    uniqueGroups = sortGroups(unique(allGroupsFinal(:)));
else
    uniqueGroups = params.groupSelection;
end

% remove measurements entries that are not in uniqueGroups
for i = 1:length(measurementsBefore)
    measurementsBefore{i} = selectMeasurements(measurementsBefore{i}, ...
        ismember(measurementsBefore{i}.groups, uniqueGroups));
end
for i = 1:length(measurementsAfter)
    measurementsAfter{i} = selectMeasurements(measurementsAfter{i}, ...
        ismember(measurementsAfter{i}.groups, uniqueGroups));
end

% use the MatrixPlotter to make a figure containing a matrix of plots
plotter = MatrixPlotter(length(uniqueGroups), params.plotterOptions{:});
plotter.axAspect = 1;
hadGuides = false(length(uniqueGroups), 1);
while plotter.next
    % find the index of the current plot
    i = plotter.index;
    
    % plot the 'before' and 'after' measurements for this plane
    hold on;
    for k = 1:length(measurementsBefore)
        % first find matching before and after measurements
        crtMeasurementsBefore = measurementsBefore{k};
        crtMeasurementsAfter = measurementsAfter{k};
        
        groupMeasurementsBefore = selectMeasurements(...
            crtMeasurementsBefore, ...
            strcmp(crtMeasurementsBefore.groups, uniqueGroups{i}));
        groupMeasurementsAfter = selectMeasurements(...
            crtMeasurementsAfter, ...
            strcmp(crtMeasurementsAfter.groups, uniqueGroups{i}));
    
        % figure mapping between indices in 'before' and 'after' measurements
        directionsBefore = groupMeasurementsBefore.directions;
        directionsAfter = groupMeasurementsAfter.directions;
        
        mappingBefore = zeros(length(directionsBefore), 1);
        mappingAfter = zeros(length(directionsAfter), 1);
        
        % XXX this is awfully inefficient, but it's unlikely to matter
        for iBef = 1:length(directionsBefore)
            crtDir = directionsBefore{iBef};
            for iAft = 1:length(directionsAfter)
                if max(abs(crtDir - directionsAfter{iAft})) < 1e-8
                    mappingBefore(iBef) = iAft;
                    mappingAfter(iAft) = iBef;
                    break;
                end
            end
        end
        
        % find directions that have a match
        hasMatchBefore = (mappingBefore ~= 0);
        hasMatchAfter = (mappingAfter ~= 0);
        
        % figure out which directions are invariant
        if isfield(crtMeasurementsAfter, 'changed')
            invariantAfter = ~groupMeasurementsAfter.changed;
        else
            invariantAfter = false(size(mappingAfter));
            invariantAfter(mappingAfter ~= 0) = ...
                abs(groupMeasurementsAfter.thresholds(hasMatchAfter) - ...
                    groupMeasurementsBefore.thresholds(mappingAfter(hasMatchAfter))) ...
                < params.invarianceTol;
        end
        invariantBefore = false(size(mappingBefore));
        invariantBefore(mappingAfter(invariantAfter)) = true;
        
        colorInvariant = selectMod(params.colorsInvariant, k);
        colorNoMatch = selectMod(params.colorsNoMatch, k);
    
        % now do the drawing
        projectedPoints = cell(1, 2);
        for j = 1:2
            switch j
                case 1
                    groupMeasurements = groupMeasurementsBefore;
                    crtColor = selectMod(params.colorsBefore, k);
                    crtMarker = selectMod(params.markersBefore, k);
                    crtSize = selectMod(params.sizesBefore, k);
                    crtColorFct = params.colorFctBefore;
                    crtHasMatch = hasMatchBefore;
                    crtInvariant = invariantBefore;
                case 2
                    groupMeasurements = groupMeasurementsAfter;
                    crtColor = selectMod(params.colorsAfter, k);
                    crtMarker = selectMod(params.markersAfter, k);
                    crtSize = selectMod(params.sizesAfter, k);
                    crtColorFct = params.colorFctAfter;
                    crtHasMatch = hasMatchAfter;
                    crtInvariant = invariantAfter;
            end

            if isfield(groupMeasurements, 'thresholdIntervals')
                errors = groupMeasurements.thresholdIntervals;
            else
                errors = [];
            end
            if ~isempty(crtColorFct)
                colors = [];
            else
                colors = repmat(crtColor(:)', ...
                    length(groupMeasurements.groups), 1);
                for cidx = 1:3
                    colors(crtInvariant, cidx) = colorInvariant(cidx);
                    colors(~crtHasMatch, cidx) = colorNoMatch(cidx);
                end
            end
            if ~isempty(groupMeasurements.groups)
                if ~hadGuides(i)
                    drawGuides = true;
                    hadGuides(i) = true;
                else
                    drawGuides = false;
                end
                projectedPoints{j} = plotTernaryThresholds(...
                    groupMeasurements, ...
                    'marker', crtMarker, ...
                    'color', colors, ...
                    'size', crtSize, ...
                    'errors', errors, ...
                    'colorFct', crtColorFct, ...
                    'drawGuides', drawGuides, ...
                    unmatched);
            end
        end
        if ~any(cellfun(@isempty, projectedPoints))
            projectedPoints{1} = projectedPoints{1}(hasMatchBefore, :);
            projectedPoints{2} = projectedPoints{2}(mappingBefore(hasMatchBefore), :);
            plot([projectedPoints{1}(:, 1) projectedPoints{2}(:, 1)]', ...
                [projectedPoints{1}(:, 2) projectedPoints{2}(:, 2)]', ...
                'color', selectMod(params.lineColors, k), 'linewidth', 1);
        end
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

%     % draw line between 'before' and 'after' predictions
%     plot([projPredBefore(:, 1) projPredAfter(:, 1)]', ...
%         [projPredBefore(:, 2) projPredAfter(:, 2)]', ...
%         'color', params.predLineColor, 'linewidth', 1);
%     
%     % overlay a neutral color over things that haven't changed at all
%     % (presumably the transformation didn't affect those directions)
%     noChangeMask = max(abs(projPredBefore - projPredAfter), [], 2) ...
%         < params.noChangeTol;
%     projNoChange = projPredBefore(noChangeMask, :);
%     plot(projNoChange(:, 1), projNoChange(:, 2), 'marker', params.predMarker, ...
%         'markersize', params.predSize', 'color', params.colorNoChange, ...
%         'linestyle', 'none');