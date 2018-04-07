function target_thresholds = matchThresholds(thresholds, groups, directions, ...
    target_groups, target_directions)
% matchThresholds Map a set of thresholds to a different ordering or
% selection of texture groups and directions.
%   target_thresholds = matchThresholds(thresholds, groups, directions, ...
%       target_groups, target_diections)
%   reorganizes the thresholds in the `thresholds` vector so that they
%   match the group names and directions from `target_groups` and
%   `target_directions`, assuming that the `thresholds` are ordered
%   according to the `groups` and `directions` cell arrays.

tol = 1e-10;

target_thresholds = nan(size(target_groups));
for i = 1:length(target_groups)
    crt_group = target_groups{i};
    crt_direction = target_directions{i};
    mask = strcmp(groups, crt_group) & ...
        cellfun(@(v) isequal(size(v), size(crt_direction)) && norm(v - crt_direction) < tol, ...
        directions);
    if sum(mask) == 0
        error([mfilename ':nomatch'], ['No match found for group ' crt_group ', direction ' num2str(crt_direction) '.']);
    end
    if sum(mask) > 1 && std(thresholds(mask)) > tol
        warning([mfilename ':multimatch'], ['Multiple distinct matches found for group ' crt_group ', direction ' num2str(crt_direction) '. Using the first one.']);
    end
    target_thresholds(i) = thresholds(find(mask, 1));
end

end