function v = ternaryextdir(group, direction)
% TERNARYEXTDIR Extend {group, direction} pair(s) to a 99-dimensional
% vector(s).
%   v = TERNARYEXTDIR(group, direction) extends the 3-component `direction`
%   in the texture group `group` to a full 99-component vector.
%
%   If `group` and `direction` are cell arrays, the function performs the
%   extension for pairs of elements, returning a cell array of 99-component
%   vectors.
%
%   This also works in the multi-group case: then the group contains two or
%   more texture group strings separated by a semicolon, and the direction
%   has 3*ngroups components.

% handle cell array input
if iscell(group) || iscell(direction)
    if ~iscell(group) || ~iscell(direction) || ~isequal(size(group) , size(direction))
        error([mfilename ':badinp'], 'The input should have the same number of groups and directions.');
    end
    v = cellfun(@(g, d) ternaryextdir(g, d), group, direction, 'uniform', false);
    return;
end

% mapping from the 99 probabilities to texture groups
coord_groups = {...
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

v = 1/3*ones(length(coord_groups), 1);

sep_groups = strtrim(strsplit(group, ';'));
for k = 1:length(sep_groups)
    crt_group = sep_groups{k};
    % get rid of a direction specification (in square brackets), if present
    sqbr_idx = find(crt_group == '[', 1);
    if ~isempty(sqbr_idx)
        crt_group = crt_group(1:sqbr_idx-1);
    end
    crt_direction = direction(3*k-2:3*k);
    crt_mask = strcmp(coord_groups, crt_group);
    v(crt_mask) = crt_direction;
end

end