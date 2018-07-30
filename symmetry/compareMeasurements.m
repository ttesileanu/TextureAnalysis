function [diff, details] = compareMeasurements(measurements1, measurements2, type, varargin)
% compareMeasurements Compare two sets of threshold measurements.
%   diff = compareMeasurements(measurements1, measurements2, 'direct')
%   returns a measure of the relative difference between the two sets of
%   measurements. The measurements structures need to contain at least
%   fields `groups`, `directions`, and `thresholds`. The calculation
%   proceeds as follows. First, the mapping between two sets is found, and
%   for every pair of matching thresholds, the difference in logarithms is
%   calculated. If a `changed` field exists in either structure, only
%   thresholds for which this field is `true` are used (though see
%   'changedOnly' option below). When the `changed` field does not exist,
%   but the 'changedOnly' option is still `true`, a heuristic is used to
%   generate the `changed` field: namely, it is assumed that any threshold
%   that changed by less than the value from the 'invarianceTol' option is
%   should have `changed == false`. The RMS of the log differences for all
%   thresholds is then calculated.
%
%   diff = compareMeasurements(measurements1, measurements2, 'group')
%   performs a similar calculation with an emphasis on the group level. Now
%   only entire groups are removed if the corresponding `changed` field is
%   `false` (so if just a few measurements inside a group are invariant,
%   they are not ignored). Also, the RMS for the log differences is first
%   calculated for every group, and then all the group values are RMSed
%   together.
%
%   [diff, details] = compareMeasurements(...) returns a structure
%   containing detailed information about the intermediate steps of the
%   calculation: the mapping between indices in the two measurement sets,
%   the log threshold differences, per-group averages, elliptic fits, etc.
%
%   Options:
%    'groupMaskFct'
%       A callable that takes in a group name and returns true to keep the
%       group and false to dismiss it (and thus not draw it).
%    'changedOnly'
%       If `true`, only measurements that have changed are considered, as
%       described above. If `false`, all measurements are used.
%    'invarianceTol'
%       Tolerance level used to define when measurements have changed in
%       cases where a `changed` field isn't present.

% XXX TODO:
%
%   diff = compareMeasurements(measurements1, measurements2, 'ellipse')
%   works only for ternary texture statistics in single-group and some
%   mixed-group planes. It fits ellipses to the measurements in each group
%   and, for those groups that aren't completely un`changed` and for which
%   fits exist in both sets of measurements, the ellipses are compared
%   according to the following metric XXX. The difference values are
%   averaged over all groups.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('groupMaskFct', @(g) true);
parser.addParameter('changedOnly', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('invarianceTol', 1e-8, @(x) isnumeric(x) && isscalar(x));

% show defaults if requested
if nargin == 1 && strcmp(measurements1, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% first find the mapping between the two sets of measurements
[details.shuffle1, details.shuffle2] = matchMeasurements(measurements1, measurements2);

% restrict to common measurements
measurements1 = selectMeasurements(measurements1, details.shuffle1 > 0);
measurements2 = selectMeasurements(measurements2, details.shuffle1(details.shuffle1 > 0));

% restrict to groups that we want to include
mask = cellfun(params.groupMaskFct, measurements1.groups);
measurements1 = selectMeasurements(measurements1, mask);
measurements2 = selectMeasurements(measurements2, mask);

details.common.measurements1 = measurements1;
details.common.measurements2 = measurements2;

% if changedOnly, find list of changed measurements
if params.changedOnly
    if isfield(measurements1, 'changed')
        changed = measurements1.changed;
    elseif isfield(measurements2, 'changed')
        changed = measurements2.changed;
    else
        changed = (abs(measurements1.thresholds - measurements2.thresholds) ...
            >= params.invarianceTol);
    end
else
    changed = true(size(measurements1.groups));
end

% ignore measurements that are NaN or infinite in either group
mask = isfinite(measurements1.thresholds) & isfinite(measurements2.thresholds);

details.common.logdiff = log(measurements1.thresholds) - log(measurements2.thresholds);

switch type
    case 'direct'
        diff = rms(details.common.logdiff(changed & mask));
    case 'group'
        if params.changedOnly
            % update 'changed' to make it consistent within groups
            uniqueGroups = sortGroups(unique(measurements1.groups));
            uniqueGroupMask = true(size(uniqueGroups));
            for i = 1:length(uniqueGroups)
                groupMask = strcmp(measurements1.groups, uniqueGroups{i});
                if any(changed(groupMask))
                    changed(groupMask) = true;
                else
                    uniqueGroupMask(i) = false;
                end
                if all(~mask(groupMask))
                    % don't include groups that have only invalid measurements
                    uniqueGroupMask(i) = false;
                end
            end
            uniqueGroups = uniqueGroups(uniqueGroupMask);
        end
        
        % calculate group averages
        details.common.uniqueGroups = uniqueGroups;
        details.common.groupDiffMeans = cellfun(@(group) ...
            rms(details.common.logdiff(strcmp(measurements1.groups, group) & mask)), ...
            uniqueGroups);
        diff = rms(details.common.groupDiffMeans);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized difference type.');
end

details.common.changed = changed;

end

function x = rms(v)
% Get root mean squared from vector.

x = sqrt(mean(v(:).^2));

end