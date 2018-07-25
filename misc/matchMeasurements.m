function [shuffle1, shuffle2] = matchMeasurements(measurements1, measurements2)
% matchMeasurements Find mapping between two sets of measurements.
%   [shuffle1, shuffle2] = matchMeasurements(measurements1, measurements2)
%   finds the mapping between the two sets of measurements, which must
%   contain fields `groups` and `directions`. Specifically, `shuffle1(i)`
%   is the index in `measurements2` corresponding to the `i`th entry in
%   `measurements1`, and conversely, `shuffle2(i)` is the index in
%   `measurements1` corresponding to the `j`th entry in `measurements2`.
%   Indices that have no match in the other measurement set are set to 0.

shuffle1 = zeros(size(measurements1.groups));
shuffle2 = zeros(size(measurements2.groups));
for i = 1:length(measurements1.groups)
    crtGroup = measurements1.groups{i};
    crtDir = measurements1.directions{i};
    
    % does this group appear in measurements2?
    groupMatchMask = strcmp(measurements2.groups, crtGroup);
    groupMatchIdxs = find(groupMatchMask);
    if ~isempty(groupMatchIdxs)
        % yes, but what about the direction?
        directionMatchSubmask = cellfun(@(d) max(abs(d - crtDir)) < 1e-6, ...
            measurements2.directions(groupMatchMask));
        nMatches = sum(directionMatchSubmask);
        if nMatches > 0
            % yes
            if nMatches > 1
                % but it appears in several locations!
                warning([mfilename ':multimatch'], ...
                    ['Entry at index %d (group %s) in measurements1 matches ' ...
                     'several entries in measurements2.'], ...
                    i, crtGroup);
            end
            submaskIdx = find(directionMatchSubmask, 1);
            matchIdx = groupMatchIdxs(submaskIdx);
            shuffle1(i) = matchIdx;
            shuffle2(matchIdx) = i;
        end
    end
end

end