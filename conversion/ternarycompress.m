function [group, direction] = ternarycompress(v)
% TERNARYCOMPRESS Compress a 99-dimensional ternary texture vector into a
% (possibly mixed) group and a direction within that group.
%   [group, direction] = TERNARYCOMPRESS(v) finds the coordinates in the
%   99-dimensional ternary texture vector `v` that are away from the origin
%   (i.e., 1/3) and identifies the group(s) that these belong to. The group
%   information is returned in `group`, with semicolons used to separate
%   the groups if several are involved. The direction within each group is
%   given by the `direction` output argument, which is a vector containing
%   3 coordinates per group, in the order in which the groups are present
%   in `group`.
%
%   For a cell array input, the function returns cell arrays of groups and
%   directions.
%
% See also: TERNARYEXTDIR.

% handle cell array input
if iscell(v)
    [group, direction] = cellfun(@(w) ternarycompress(w), v, 'uniform', false);
    return;
end

% mapping from the 99 probabilities to texture groups
coordGroups = {...
    'A_1'       'A_1'       'A_1'       'AB_1_1'    'AB_1_1'    'AB_1_1' ...
    'AB_1_2'    'AB_1_2'    'AB_1_2'    'AC_1_1'    'AC_1_1'    'AC_1_1'  ...
    'AC_1_2'    'AC_1_2'    'AC_1_2'    'BC_1_1'    'BC_1_1'    'BC_1_1'  ...
    'BC_1_2'    'BC_1_2'    'BC_1_2'    'AD_1_1'    'AD_1_1'    'AD_1_1'  ...
    'AD_1_2'    'AD_1_2'    'AD_1_2'    'ABC_1_1_1' 'ABC_1_1_1' 'ABC_1_1_1' ...
    'ABC_1_2_2' 'ABC_1_2_2' 'ABC_1_2_2' 'ABC_1_2_1' 'ABC_1_2_1' 'ABC_1_2_1' ...
    'ABC_1_1_2' 'ABC_1_1_2' 'ABC_1_1_2' 'ABD_1_1_1' 'ABD_1_1_1' 'ABD_1_1_1' ...
    'ABD_1_2_2' 'ABD_1_2_2' 'ABD_1_2_2' 'ABD_1_2_1' 'ABD_1_2_1' 'ABD_1_2_1' ...
    'ABD_1_1_2' 'ABD_1_1_2' 'ABD_1_1_2' 'ACD_1_1_1' 'ACD_1_1_1' 'ACD_1_1_1' ...
    'ACD_1_2_2' 'ACD_1_2_2' 'ACD_1_2_2' 'ACD_1_2_1' 'ACD_1_2_1' 'ACD_1_2_1' ...
    'ACD_1_1_2' 'ACD_1_1_2' 'ACD_1_1_2' 'BCD_1_1_1' 'BCD_1_1_1' 'BCD_1_1_1' ...
    'BCD_1_2_2' 'BCD_1_2_2' 'BCD_1_2_2' 'BCD_1_2_1' 'BCD_1_2_1' 'BCD_1_2_1' ...
    'BCD_1_1_2' 'BCD_1_1_2' 'BCD_1_1_2' ...
    'ABCD_1_1_1_1'    'ABCD_1_1_1_1'    'ABCD_1_1_1_1' ...
    'ABCD_1_2_2_2'    'ABCD_1_2_2_2'    'ABCD_1_2_2_2' ...
    'ABCD_1_2_1_1'    'ABCD_1_2_1_1'    'ABCD_1_2_1_1' ...
    'ABCD_1_1_2_2'    'ABCD_1_1_2_2'    'ABCD_1_1_2_2' ...
    'ABCD_1_1_2_1'    'ABCD_1_1_2_1'    'ABCD_1_1_2_1' ...
    'ABCD_1_2_1_2'    'ABCD_1_2_1_2'    'ABCD_1_2_1_2' ...
    'ABCD_1_2_2_1'    'ABCD_1_2_2_1'    'ABCD_1_2_2_1' ...
    'ABCD_1_1_1_2'    'ABCD_1_1_1_2'    'ABCD_1_1_1_2'};

% find entries in v that are not at the origin, within some tolerance
tol = 1e-6;
changed = abs(v - 1/3) > tol;

% find the groups
groups = unique(coordGroups(changed), 'stable');

% find the direction information in changed groups
groupMask = ismember(coordGroups, groups);
direction = v(groupMask);

% if it's a pair of two groups, check whether it's within a nice 2-d
% subspace
if length(groups) == 2
    try
        [~, groupDirs] = ternary6tomix2(direction);
    catch
        groupDirs = [];
    end
    if ~isempty(groupDirs)
        groups = arrayfun(@(i) [groups{i} '[' int2str(groupDirs(i)) ']'], ...
            1:length(groups), 'uniform', false);
    end
end

% collapse the groups to a single string
group = cell2mat(flatten([groups ; [repmat({';'}, 1, length(groups)-1) {''}]])');

end