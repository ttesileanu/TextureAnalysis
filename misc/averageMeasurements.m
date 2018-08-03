function measurements = averageMeasurements(measurements1, measurements2, varargin)
% averageMeasurements Calculate average measurements from two measurement
% sets.
%   measurements = averageMeasurements(measurements1, measurements2)
%   calculates the average `thresholds` and `thresholdIntervals` for
%   measurements in the two measurement sets.
%
%   Options:
%    'skipMissing'
%       If `true`, thresholds that exist in one structure but not in the
%       other are skipped. If `false`, all thresholds that exist in at
%       least one structure are kept.
%    'meanFct'
%       Function to use for averaging. This should have the same interface
%       as Matlab's `mean` function.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('skipMissing', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('meanFct', @geomean);

% show defaults if requested
if nargin == 1 && strcmp(measurements1, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% find matching measurements
[shuffle1, shuffle2] = matchMeasurements(measurements1, measurements2);
mask1 = (shuffle1 ~= 0);
mask2 = (shuffle2 ~= 0);

% start creating output structure by using the common measurements
measurements.groups = measurements1.groups(mask1);
measurements.directions = measurements1.directions(mask1);
if isfield(measurements1, 'nSubjects')
    measurements.nSubjects = measurements1.nSubjects(mask1);
else
    measurements.nSubjects = zeros(sum(mask1), 1);
end
if isfield(measurements2, 'nSubjects')
    measurements.nSubjects = measurements.nSubjects ...
        + measurements2.nSubjects(shuffle1(mask1));
end
measurements.multi = (cellfun(@(g) sum(g == ';'), measurements.groups) > 0);
measurements.thresholds = params.meanFct(...
    [measurements1.thresholds(mask1) measurements2.thresholds(shuffle1(mask1))], 2);
if isfield(measurements1, 'thresholdIntervals') && isfield(measurements2, 'thresholdIntervals')
    lo = params.meanFct([measurements1.thresholdIntervals(mask1, 1), ...
        measurements2.thresholdIntervals((shuffle1(mask1)), 1)], 2);
    hi = params.meanFct([measurements1.thresholdIntervals(mask1, 2), ...
        measurements2.thresholdIntervals((shuffle1(mask1)), 2)], 2);
    measurements.thresholdIntervals = [lo hi];
elseif isfield(measurements1, 'thresholdIntervals')
    measurements.thresholdIntervals = measurements1.thresholdIntervals(mask1);
elseif isfield(measurements2, 'thresholdIntervals')
    measurements.thresholdIntervals = measurements2.thresholdIntervals((shuffle1(mask1)));
end

% now concatenate measurements that only appear in set 1
if sum(~mask1) > 0
    measurements = catMeasurements(measurements, ...
        selectMeasurements(measurements1, ~mask1));
end

% finally concatenate measurements that only appear in set 2
if sum(~mask2) > 0
    measurements = catMeasurements(measurements, ...
        selectMeasurements(measurements2, ~mask2));
end

end
