function shuffled = shuffleMeasurements(measurements, varargin)
% shuffleMeasurements Randomly shuffle thresholds.
%   shuffled = shuffleMeasurements(measurements) randomly shuffles the
%   thresholds in the `measurements` structure. If a `thresholdIntervals`
%   field is present, it is shuffled in the same way.
%
%   Options:
%    'group'
%       If `true`, shuffle thresholds only within groups.
%    'cyclic'
%       If `true`, generate only cyclic permutation (either overall, or
%       within each group, if 'groupShuffle' is `true`).

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('group', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('cyclic', false, @(b) islogical(b) && isscalar(b));

% show defaults if requested
if nargin == 1 && strcmp(measurements, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

if ~params.group
    if ~params.cyclic
        permutation = randperm(length(measurements.thresholds));
    else
        permutation = circshift(1:length(measurements.thresholds), ...
            randi(length(measurements.thresholds)));
    end
else
    uniqueGroups = unique(measurements.groups);
    permutation = 1:length(measurements.thresholds);
    for i = 1:length(uniqueGroups)
        mask = strcmp(measurements.groups, uniqueGroups{i});
        crtN = sum(mask);
        if ~params.cyclic
            crtPerm = randperm(crtN);
        else
            crtPerm = circshift(1:crtN, randi(crtN));
        end
        crtIdxs = permutation(mask);
        permutation(mask) = crtIdxs(crtPerm);
    end
end

shuffled = measurements;
shuffled.thresholds = shuffled.thresholds(permutation);
if isfield(shuffled, 'thresholdIntervals')
    shuffled.thresholdIntervals = shuffled.thresholdIntervals(permutation, :);
end

end
