function [projectedThresholds] = plotTernaryThresholds(...
    thresholds, varargin)
% plotTernaryThresholds Plot thresholds in a single- or double-group
% ternary texture plane.
%   plotTernaryThresholds(thresholds, directions) plots the given
%   `thresholds` in a ternary texture plane, assuming that each threshold
%   corresponds to a direction given in the `directions` cell array.
%   Whether the plane should be single- or double-group is inferred from
%   the size of the elements of the `directions` array. Error bars can also
%   be drawn using the 'errors' option (see below).
%
%   plotTernaryThresholds(measurements) uses a structure with fields
%   `thresholds` and `directions` to make the plot. A `groups` field, if it
%   exists, can be used to label the axes.
%
%   projectedThresholds = plotTernaryThresholds(...) returns the calculated
%   2d positions of the threholds in the texture plane.
%
%   Options:
%    'marker'
%       Marker to use.
%    'color'
%       3-component RGB vector giving the color of the symbols to use. This
%       can also be a matrix in which every row gives the color of a
%       matching entry in `thresholds` and `directions`.
%    'size'
%       Symbol size to use. This can also be a vector in which every entry
%       gives the size for a matching entry in `thresholds` and
%       `directions`.
%    'errors'
%       Matrix of error intervals for the thresholds, in a format where
%       each row reprents [lo, hi] limits.
%    'errorColor'
%       3-component RGB vector giving the color with which to mix the symbol
%       color to obtain the error bar color.
%    'errorBars'
%       Set to true to draw error bars, false otherwise.
%    'ellipse'
%       Set to true to show a best-fit ellipse, false otherwise.
%    'colorFct'
%       Callable giving the color to use for each subject. This should
%       behave like a function, `colorFct(subjectName)` yielding an RGB
%       color. If given, this overrides the 'color' option.
%    'drawGuides'
%       Set to `true`, `false`, or 'auto'. If `true`, or if set to 'auto'
%       and `ishold` is false, the function draws appropriate guide lines
%       (using either `drawTernaryTriangle` or `drawTernaryMixedBackground`)
%       before plotting the thresholds.
%    'triangleOptions'
%       Options to pass to drawTernaryTriangle.
%    
%   See also: plotTernaryMatrix.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('directions', {}, @(c) iscell(c));
parser.addParameter('marker', 'x', @(s) ischar(s));
parser.addParameter('color', [0.8 0.3 0.3], @(v) isempty(v) || ...
    (ismatrix(v) && isnumeric(v) && size(v, 2) == 3));
parser.addParameter('size', 5, @(x) isnumeric(x) && isvector(x) && all(x > 0));
parser.addParameter('errors', [], @(v) isempty(v) || (ismatrix(v) && size(v, 2) == 2 && isnumeric(v)));
parser.addParameter('errorColor', [0.3 0.3 0.3], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('errorBars', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('ellipse', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('colorFct', []);
parser.addParameter('drawGuides', 'auto', @(b) (islogical(b) && isscalar(b)) || strcmp(b, 'auto'));
parser.addParameter('triangleOptions', {}, @(c) iscell(c) && isvector(c));

% show defaults if requested
if nargin == 1 && strcmp(thresholds, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% handle errorBars
if ~params.errorBars
    params.errors = [];
end

% handle struct input
if isstruct(thresholds)
    measurements = thresholds;
    thresholds = measurements.thresholds;
    directions = measurements.directions;
else
    measurements = struct;
    directions = params.directions;
end

% handle default drawGuides
if strcmp(params.drawGuides, 'auto')
    params.drawGuides = ~ishold;
end

% don't display empty inputs
if isempty(thresholds)
    projectedThresholds = [];
    return;
end

% figure out how whether this is a single-group or a mixed-group plane
if mod(length(directions{1}), 3) ~= 0
    error([mfilename ':baddirsz'], 'Directions vector sizes aren''t compatible with ternary planes.');
end
nGroups = length(directions{1}) / 3;
if any(cellfun(@(c) length(c) ~= nGroups*3, directions))
    error([mfilename ':mixeddirs'], 'Not all directions have the same size.');
end

% do we need to draw guides?
if params.drawGuides
    if isfield(measurements, 'groups') && ~isempty(measurements.groups)
        setupTernaryGuide(measurements.groups{1}, 'triangleOptions', params.triangleOptions);
    else
        setupTernaryGuide(nGroups, 'triangleOptions', params.triangleOptions);
    end
end

% convert predictions, thresholds, and threshold intervals to the 2d
% coordinates necessary for plotting
switch nGroups
    case 1
        % single-group plane
        projectedThresholds = ternary3to2(ternaryrec(thresholds, directions));
        if ~isempty(params.errors)
            projectedLo = ternary3to2(ternaryrec(params.errors(:, 1), directions));
            projectedHi = ternary3to2(ternaryrec(params.errors(:, 2), directions));
        end
    case 2
        % mixed-group plane
        projectedThresholds = ternary6tomix2(ternaryrec(thresholds, directions));
        if ~isempty(params.errors)
            projectedLo = ternary6tomix2(ternaryrec(params.errors(:, 1), directions));
            projectedHi = ternary6tomix2(ternaryrec(params.errors(:, 2), directions));
        end
    otherwise
        error([mfilename ':badngrp'], 'This function only works with single groups and pairs of groups.');
end

% figure out colors
if isempty(params.colorFct)
    colors = params.color;
else
    if isfield(measurements, 'subjects')
        colors = cell2mat(cellfun(@(s) flatten(params.colorFct(s))', ...
            measurements.subjects(:), 'uniform', false));
    else
        colors = cell2mat(arrayfun(@(th) flatten(params.colorFct('subject'))', ...
            thresholds(:), 'uniform', false));
    end
end

% draw error bars for measured data, if we have them
if ~isempty(params.errors)
    % don't plot errorbars that are not finite, though
    errMask = all(isfinite(projectedLo), 2) & all(isfinite(projectedHi), 2);
    
    errorColor = mixcolor(colors, params.errorColor);
    
    projLoMasked = projectedLo(errMask, :);
    projHiMasked = projectedHi(errMask, :);
    
    if ~isempty(projLoMasked)        
        if size(errorColor, 1) > 1
            errColMasked = errorColor(errMask, :);
        else
            errColMasked = errorColor;
        end
        
        % argh: need to draw each line separately if colors change
        if size(errColMasked, 1) == 1
            plot([projLoMasked(:, 1)' ; projHiMasked(:, 1)'], ...
                [projLoMasked(:, 2)' ; projHiMasked(:, 2)'], ...
                'color', errColMasked, 'linewidth', 0.5);
        else
            for k = 1:size(errColMasked, 1)
                plot([projLoMasked(k, 1) projHiMasked(k, 1)], ...
                    [projLoMasked(k, 2) projHiMasked(k, 2)], ...
                    'color', errColMasked(k, :), 'linewidth', 0.5);
            end
        end
    end
end

% draw measured data
h = scatter(projectedThresholds(:, 1), projectedThresholds(:, 2), ...
    params.size.^2, colors, params.marker(1));
h.LineWidth = 1;

% draw ellipses
if params.ellipse
    % measurements
    finiteMask = all(isfinite(projectedThresholds), 2);
    if sum(finiteMask) > 2
        % remove non-finite entries
        measM = fitEllipse(projectedThresholds(finiteMask, :));
        avgColor = mean(colors, 1);
        ellipse(0, 0, measM, 'color', mixcolor(avgColor, [0.5 0.5 0.5]));
    end
end

end
