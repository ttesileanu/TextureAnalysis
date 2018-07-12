function plotTernaryThresholds(thresholds, directions, varargin)
% plotTernaryThresholds Plot thresholds in a single- or double-group
% ternary texture plane.
%   plotTernaryThresholds(thresholds, directions) plots the given
%   `thresholds` in a ternary texture plane, assuming that each threshold
%   corresponds to a direction given in the `directions` cell array.
%   Whether the plane should be single- or double-group is inferred from
%   the size of the elements of the `directions` array. Error bars can also
%   be drawn using the 'errors' option (see below).
%
%   Options:
%    'marker'
%       Marker to use in general.
%    'markerNoError'
%       marker to use if error bars are provided, but they are infinite or
%       NaN for a particular measurement.
%    'color'
%       3-component RGB vector giving the color of the symbols to use.
%    'size'
%       Symbol size to use.
%    'errors'
%       Matrix of error intervals for the thresholds, in a format where
%       each row reprents [lo, hi] limits.
%    'errorColor'
%       3-component RGB vector giving the color of the error bars.
%    'ellipse'
%       Set to true to show a best-fit ellipse, false otherwise.
%    
%   See also: plotTernaryPlane, plotTernaryMatrix.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('marker', 'x', @(s) ischar(s));
parser.addParameter('markerNoError', 'o', @(s) ischar(s));
parser.addParameter('color', [0.8 0.3 0.3], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('size', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('errors', [], @(v) isempty(v) || (ismatrix(v) && size(v, 2) == 2 && isnumeric(v)));
parser.addParameter('errorColor', [0.5 0.5 0.5], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('ellipse', false, @(b) islogical(b) && isscalar(b));

% show defaults if requested
if nargin == 1 && strcmp(thresholds, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% figure out how whether this is a single-group or a mixed-group plane
if mod(directions{1}, 3) ~= 0
    error([mfilename ':baddirsz'], 'Directions vector sizes aren''t compatible with ternary planes.');
end
nGroups = length(directions{1}) / 3;
if any(cellfun(@(c) length(c) ~= nGroups*3, directions))
    error([mfilename ':mixeddirs'], 'Not all directions have the same size.');
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

% draw error bars for measured data, if we have them
if ~isempty(params.errors)
    % don't plot errorbars that are not finite, though
    errMask = all(isfinite(projectedLo), 2) & all(isfinite(projectedHi), 2);
    
    plot([projectedLo(errMask, 1)' ; projectedHi(errMask, 1)'], ...
         [projectedLo(errMask, 2)' ; projectedHi(errMask, 2)'], ...
         'color', params.errorColor, 'linewidth', 0.5);
end

% draw measured data
if ~isempty(params.errors)
    % draw points without error bars as small open circles
    circleMask = ~errMask;
else
    circleMask = false(size(thresholds));
end

% draw points with error bars
plot(projectedThresholds(~circleMask, 1), ...
    projectedThresholds(~circleMask, 2), params.marker, 'linewidth', 1, ...
    'markersize', params.size, 'color', params.color);
% draw points without error bars
plot(projectedThresholds(circleMask, 1), ...
    projectedThresholds(circleMask, 2), params.markerNoError, 'linewidth', 1, ...
    'markersize', params.size/2, 'color', params.color);

% draw ellipses
if params.ellipse
    % measurements
    finiteMask = all(isfinite(projectedThresholds), 2);
    if sum(finiteMask) > 2
        % remove non-finite entries
        measM = fitEllipse(projectedThresholds(finiteMask, :));
        ellipse(0, 0, inv(measM), 'color', mixcolor(params.color, [0.5 0.5 0.5]));
    end
end

end