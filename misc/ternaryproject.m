function w = ternaryproject(v, group)
% TERNARYPROJECT Project a 99-dimensional ternary texture vector onto the
% directions corresponding to a given group or set of groups.
%   w = TERNARYPROJECT(v, group) returns a 3d vector that represents the
%   projection of the 99-dimensional ternary texture vector `v` onto the
%   given texture group. If `v` is a matrix, each column is transformed
%   separately, returning a matrix `w`.
%
%   w = TERNARYPROJECT(v, {group1, ..., groupk}) projects onto the
%   3k-dimensional space spanned by the `k` groups.
%
%   Alternatively, a multi-group projection can be obtained by separating
%   the different groups by semicolon, like so: 'AB_1_1 ; AC_1_2'. Separate
%   directions of groups can also be used, so that 'AB_1_2[1]' refers to
%   the second (mid-level gray) direction in the AB_1_2 plane.
%
%   See also: TERNARYEXTDIR.

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

% handle vector vs. matrix
if isvector(v)
    v = v(:);
end
nv = size(v, 2);

% handle single groups, semicolon-separated multi-groups, and cell array
% multi-groups
if ischar(group)
    sep_groups = strtrim(strsplit(group, ';'));
else
    sep_groups = group;
end

k = length(sep_groups);
w = zeros(3*k, nv);

crt = 1;
for k = 1:length(sep_groups)
    crt_group = sep_groups{k};
    % handle direction specification (in square brackets), if present
    sqbr_idx = find(crt_group == '[', 1);
    if ~isempty(sqbr_idx)
        crt_dir = str2double(crt_group(sqbr_idx+1));
        crt_group = crt_group(1:sqbr_idx-1);
        crt_idxs = find(strcmp(coord_groups, crt_group));
        crt_idx = crt_idxs(1 + crt_dir);
        w(crt, :) = v(crt_idx, :);
        crt = crt + 1;
    else
        crt_mask = strcmp(coord_groups, crt_group);
        w(crt:crt+2, :) = v(crt_mask, :);
        crt = crt + 3;
    end
end
w = w(1:crt-1, :);

end