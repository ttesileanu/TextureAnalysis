function group = normalizegroup(group)
% NORMALIZEGROUP Normalize texture group name.
%   group = NORMALIZEGROUP(group) normalizes the group name by, e.g.,
%   padding it with spaces so that it sorts in a natural way and can be
%   used for comparisons.
%
%   If `group` is a cell array, the function normalizes each element. The
%   function also works for multi-group names (where the groups are
%   separated by semicolons).

% handle empty input
if isempty(group)
    return;
end

% handle cell array version
if iscell(group)
    group = cellfun(@normalizegroup, group, 'uniform', false);
    return;
end

% handle semicolon-separated multi-group names
multi_group = strtrim(strsplit(group, ';'));
if length(multi_group) > 1
    multi_group = normalizegroup(multi_group);
    parts = cell(1, 2*length(multi_group)-1);
    parts(1:2:end) = multi_group;
    parts(2:2:end) = {';'};
    group = strcat(parts{:});
    return;
end

% remove white space
group = group(~isspace(group));

under_start = find(group == '_', 1);
letters0 = group(1:under_start-1);
letters = pad(letters0, 4, 'left');

dir_start = find(group == '[', 1);
if ~isempty(dir_start)
    dir_end = find(group == ']', 1);
    direction = group(dir_start:dir_end);
else
    dir_start = length(group) + 1;
    direction = '[0]';
end

under0 = group(under_start:dir_start-1);
under = [repmat('_0', 1, 4-length(letters0)) under0];

group = [letters under direction];

end