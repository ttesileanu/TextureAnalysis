function [difference, details] = compareMeasurements(measurements1, measurements2, type, varargin)
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
%   diff = compareMeasurements(measurements1, measurements2, 'ellipse')
%   works only for ternary texture statistics in single-group and some
%   mixed-group planes. It fits ellipses to the measurements in each group
%   and, for those groups that aren't completely un`changed` and for which
%   fits exist in both sets of measurements, the ellipses are compared
%   according to the following metric. Each ellipse is first converted to a
%   3d vector
%       wvec = (a+b)/2 * (ecc*e1 + sqrt(1-ecc^2)*zhat) ,
%   where a and b are the lengths of the two semiaxes, e1 is the unit
%   vector pointing in the direction of the long semiaxes, zhat is the unit
%   vector in the z direction, and ecc is the eccentricity,
%       ecc = sqrt(1 - b^2/a^2) .
%   The difference between the two sets of measurements is the difference
%   between their associated 3d vectors. The difference values for all
%   planes where they could be calculated are RMSed together.
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
%    'normalize'
%       Before comparison, normalize both sets of measurements by
%       subtracting out, from each set, the median of the log threshold.
%       This only applies to the 'direct' and 'group' methods.
%    'hiLoRatioLimit'
%       For measurements that have error bars (in the form of a field
%       called `thresholdIntervals`), only thresholds for which the ratio
%       between the 'hi' and 'lo' limits is below the value provided here
%       are considered for the calculation. The ratios from the first
%       measurements structure are used.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('groupMaskFct', @(g) true);
parser.addParameter('changedOnly', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('invarianceTol', 1e-8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('normalize', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('hiLoRatioLimit', inf, @(x) isnumeric(x) && isscalar(x));

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

if isfield(measurements1, 'thresholdIntervals') && isfinite(params.hiLoRatioLimit)
    % also ignore measurements that have too high hi-to-lo ratio
    logRatios = diff(log(measurements1.thresholdIntervals), [], 2);
    mask = mask & (logRatios <= log(params.hiLoRatioLimit));
end

logThresholds1 = log(measurements1.thresholds);
logThresholds2 = log(measurements2.thresholds);
if params.normalize
    logThresholds1 = logThresholds1 - nanmedian(logThresholds1);
    logThresholds2 = logThresholds2 - nanmedian(logThresholds2);
end
details.common.logdiff = logThresholds2 - logThresholds1;

switch type
    case 'direct'
        difference = rms(details.common.logdiff(changed & mask));
        details.common.nAveraged = sum(changed & mask);
    case 'group'
        uniqueGroups = sortGroups(unique(measurements1.groups));
        uniqueGroupMask = true(size(uniqueGroups));
        for i = 1:length(uniqueGroups)
            groupMask = strcmp(measurements1.groups, uniqueGroups{i});
            if params.changedOnly
                if any(changed(groupMask))
                    % update 'changed' to make it consistent within groups
                    changed(groupMask) = true;
                else
                    uniqueGroupMask(i) = false;
                end
            end
            if all(~mask(groupMask))
                % don't include groups that have only invalid measurements
                uniqueGroupMask(i) = false;
            end
        end
        uniqueGroups = uniqueGroups(uniqueGroupMask);
        
        % calculate group averages
        details.common.uniqueGroups = uniqueGroups;
        details.common.groupDiffMeans = cellfun(@(group) ...
            rms(details.common.logdiff(strcmp(measurements1.groups, group) & mask)), ...
            uniqueGroups);
        difference = rms(details.common.groupDiffMeans);
    case 'ellipse'
        uniqueGroups = sortGroups(unique(measurements1.groups));
        uniqueGroupMask = true(size(uniqueGroups));
        wvecs1 = cell(size(uniqueGroups));
        wvecs2 = cell(size(uniqueGroups));
        wvecDiff = nan(size(uniqueGroups));
        for i = 1:length(uniqueGroups)
            groupMask = strcmp(measurements1.groups, uniqueGroups{i});
            if params.changedOnly
                if any(changed(groupMask))
                    % update 'changed' to make it consistent within groups
                    changed(groupMask) = true;
                else
                    uniqueGroupMask(i) = false;
                    continue;
                end
            end
            if all(~mask(groupMask))
                % don't include groups that have only invalid measurements
                uniqueGroupMask(i) = false;
                continue;
            end
            
            % project to a texture plane -- make sure these are ternary stats
            nGroups = 1 + sum(uniqueGroups{i} == ';');
            crtDirections = measurements1.directions(groupMask);
            nGray = length(crtDirections{1}) / nGroups;
            if nGray ~= 3
                error([mfilename ':badng'], 'Ellipse mode only works with ternary statistics.');
            end
            
            crtThresholds1 = measurements1.thresholds(groupMask);
            crtThresholds2 = measurements2.thresholds(groupMask);
            
            crtThreshVectors1 = ternaryrec(crtThresholds1, crtDirections);
            crtThreshVectors2 = ternaryrec(crtThresholds2, crtDirections);
            
            if nGroups == 1
                projectionFct = @ternary3to2;
            elseif nGroups == 2
                projectionFct = @ternary6tomix2;
            else
                error([mfilename ':badord'], 'Ellipse mode only works with single groups or pairs of groups.');
            end
            
            crtProjected1 = projectionFct(crtThreshVectors1);
            crtProjected2 = projectionFct(crtThreshVectors2);
            
            % fit ellipses
            finiteMask1 = all(isfinite(crtProjected1), 2);
            finiteMask2 = all(isfinite(crtProjected2), 2);
            M1 = fitEllipse(crtProjected1(finiteMask1, :));
            M2 = fitEllipse(crtProjected2(finiteMask2, :));
            
            % calculate and store 3d vectors
            wvecs1{i} = ellipseTo3d(M1);
            wvecs2{i} = ellipseTo3d(M2);
            
            % don't compare if either ellipse wasn't fit
            if ~all(isfinite(wvecs1{i})) || ~all(isfinite(wvecs2{i}))
                uniqueGroupMask(i) = false;
                continue;
            end
            
            % calculate difference
            wvecDiff(i) = norm(wvecs1{i} - wvecs2{i});
        end
        
        % keep only the groups for which we managed to calculate a difference
        uniqueGroups = uniqueGroups(uniqueGroupMask);
        wvecs1 = wvecs1(uniqueGroupMask);
        wvecs2 = wvecs2(uniqueGroupMask);
        wvecDiff = wvecDiff(uniqueGroupMask);
        
        % store the w vectors
        details.common.uniqueGroups = uniqueGroups;
        details.common.ellipseVectors1 = wvecs1;
        details.common.ellipseVectors2 = wvecs2;
        details.common.ellipseDifferences = wvecDiff;
        
        % calculate overall difference
        difference = rms(details.common.ellipseDifferences);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized difference type.');
end

details.common.changed = changed;

end

function x = rms(v)
% Get root mean squared from vector.

x = sqrt(mean(v(:).^2));

end

function wvec = ellipseTo3d(M)
% Calculate 3d vector corresponding to an ellipse, when the ellipse is
% identified by its matrix M such that x'*M*x = 1.

if any(~isfinite(M(:))) || any(eig(M) < -eps)
    wvec = nan(3, 1);
    return;
end

% find eccentricity, average radius, and direction of large semiaxis
[V, D] = eigsorted(M);
a = 1./sqrt(D(1, 1));
b = 1./sqrt(D(2, 2));
avgRad = (a + b) / 2;
ecc = sqrt(1 - b^2 / a^2);
e1 = real(V(:, 1));

% extend to 3d
zhat = [0 0 1]';
e1ext = [e1 ; 0];

wvec = avgRad*(ecc*e1ext + sqrt(1 - ecc^2)*zhat);

end
