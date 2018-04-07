function sorted = sortgroups(groups)
% SORTGROUPS Sort a cell array of group names in a natural order.
%   sorted = SORTGROUPS(groups) takes in a cell array of group names and
%   sorts them in a pleasing way. Lower-order groups are sorted before
%   higher-order ones, then sorting is done by the letter part of the group
%   (e.g., AB comes before AC), then by the multipliers (e.g., ABCD_1_1_2_1
%   comes before ABCD_1_2_2_1). For multi-group names (in which group
%   strings are separated by a semicolon ';'), sorting is done first by the
%   first group, then by the second group, etc. When groups contain
%   direction indications, they are the last key to be sorted by (e.g.,
%   AC_1_1[0] comes before AC_1_1[1], but after AB_1_1[2]). Single groups
%   come before pairs of groups, which come before triples, etc.

% to perform the sorting, we normalize the naming of the groups
sorted = groups;

% first remove all white space
groups = cellfun(@(s) s(~isspace(s)), groups, 'uniform', false);

group_count = zeros(length(groups), 1);
for i = 1:length(groups)
    crt_groups = strsplit(groups{i}, ';');
    group_count(i) = length(crt_groups);
    crt_normalized = cellfun(@normalize, crt_groups, 'uniform', false);
    crt_norm_with_sep = cell(1, 2*group_count(i)-1);
    crt_norm_with_sep(1:2:end) = crt_normalized;
    crt_norm_with_sep(2:2:end) = {';'};
    groups{i} = cell2mat(crt_norm_with_sep);
end

[~, order] = sort(groups);
sorted = sorted(order);
group_count = group_count(order);

% now make sure we put smaller tuples of groups before larger ones
[~, order2] = sort(group_count);
sorted = sorted(order2);

end

function group = normalize(group)
% Normalize a group name.

under_start = find(group == '_', 1);
letters0 = group(1:under_start-1);
letters = pad(letters0, 4, 'left');

dir_start = find(group == '[', 1);
if ~isempty(dir_start)
    dir_end = find(group == ']', 1);
    direction = group(dir_start:dir_end);
else
    dir_start = length(group) + 1;
    direction = '[ ]';
end

under0 = group(under_start:dir_start-1);
under = [repmat('_0', 1, 4-length(letters0)) under0];

group = [letters under direction];

end