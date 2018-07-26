function [sorted, shuffle] = sortGroups(groups)
% sortGroups Sort a cell array of group names in a natural order.
%   sorted = sortGroups(groups) takes in a cell array of group names and
%   sorts them in a pleasing way. Lower-order groups are sorted before
%   higher-order ones, then sorting is done by the letter part of the group
%   (e.g., AB comes before AC), then by the multipliers (e.g., ABCD_1_1_2_1
%   comes before ABCD_1_2_2_1). For multi-group names (in which group
%   strings are separated by a semicolon ';'), sorting is done first by the
%   first group, then by the second group, etc. When groups contain
%   direction indications, they are the last key to be sorted by (e.g.,
%   AC_1_1[0] comes before AC_1_1[1], but after AB_1_1[2]). Single groups
%   come before pairs of groups, which come before triples, etc.
%
%   [sorted, shuffle] = sortGroups(groups) also returns the ordering, such
%   that `sorted = groups(shuffle)`.

% to perform the sorting, we normalize the naming of the groups
sorted = groups;
groups = normalizeGroup(groups);

% find multi-groups
groupCount = cellfun(@(s) sum(s == ';'), groups);

[~, order] = sort(groups);
sorted = sorted(order);
groupCount = groupCount(order);

% now make sure we put smaller tuples of groups before larger ones
[~, order2] = sort(groupCount);
sorted = sorted(order2);

shuffle = order(order2);

end