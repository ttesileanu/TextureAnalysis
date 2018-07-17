function groups = getCoordinateMapping(G, type)
% getCoordinateMapping Find the mapping from texture coordinates to group
% names.
%   groups = getCoordinateMapping(G) returns the canonical order of texture
%   groups with `G` gray levels. For `G >= 3`, this matches the order in
%   which the groups appear in the parametrization used by `analyzeTexture`
%   (but each group is listed just once). This is not true for `G == 2`.
%   For all `G`, the group ordering returned here matches the
%   parametrization obtained from `expandTextureStats`. See more options
%   below.
%
%   getCoordinateMapping(G, 'analyzeFull') returns a list with the same
%   size as the output from `analyzeTexture`. This means that each group
%   appears `G-1` times. This also respects the `analyzeTexture`
%   parametrization even for `G == 2`.
%
%   getCoordinateMapping(G, 'analyzeShort') uses the `analyzeTexture`
%   ordering, but lists every group only once.
%
%   getCoordinateMapping(G, 'expanded') returns a list with the same size
%   as the output from `expandTextureStats`, that is, each group appears
%   `G` times. The ordering is that from `expandTextureStats` also.
%
%   getCoordinateMapping(G, 'unique') is the default.
%
%   This function uses Jonathan's framework for G >= 3.
%
% See also: analyzeTexture.

% set default type
if nargin < 2
    type = 'unique';
end

% first figure out the names of the unique groups
if G == 2
    % G == 2 needs special treatment for historical reasons
    switch type
        case {'unique', 'expanded'}
            uniqueGroups = {'A_1'    'AB_1_1'    'AC_1_1'    'BC_1_1'    'AD_1_1' ...
                'ABC_1_1_1'    'ABD_1_1_1'    'ACD_1_1_1'    'BCD_1_1_1'    'ABCD_1_1_1_1'};
        case {'analyzeFull', 'analyzeShort'}
            uniqueGroups = {'A_1', 'AC_1_1', 'AB_1_1', 'AD_1_1', 'BC_1_1', ...
                'ABC_1_1_1', 'BCD_1_1_1', 'ABD_1_1_1', 'ACD_1_1_1', 'ABCD_1_1_1_1'};
        otherwise
            error([mfilename ':badtype'], 'Unrecognized choice for return type.');
    end
else
    % hard-coding G == 3 to reduce dependence on Jonathan's framework
    if G == 3
        uniqueGroups = {...
            'A_1'           'AB_1_1'        'AB_1_2'        'AC_1_1'  ...
            'AC_1_2'        'BC_1_1'        'BC_1_2'        'AD_1_1'  ...
            'AD_1_2'        'ABC_1_1_1'     'ABC_1_2_2'     'ABC_1_2_1' ...
            'ABC_1_1_2'     'ABD_1_1_1'     'ABD_1_2_2'     'ABD_1_2_1' ...
            'ABD_1_1_2'     'ACD_1_1_1'     'ACD_1_2_2'     'ACD_1_2_1' ...
            'ACD_1_1_2'     'BCD_1_1_1'     'BCD_1_2_2'     'BCD_1_2_1' ...
            'BCD_1_1_2'     'ABCD_1_1_1_1'  'ABCD_1_2_2_2'  'ABCD_1_2_1_1' ...
            'ABCD_1_1_2_2'  'ABCD_1_1_2_1'  'ABCD_1_2_1_2'  'ABCD_1_2_2_1' ...
            'ABCD_1_1_1_2'};
    else
        mtc = analyzeTexture('mtc', G);
        uniqueGroups = cellfun(@(s) s.name, mtc.coord_groups(:)', 'uniform', false);
    end
end

% now repeat groups if necessary
switch type
    case 'analyzeFull'
        nReps = G-1;
    case 'expanded'
        nReps = G;
    otherwise
        nReps = 1;
end

if nReps > 1
    groups = flatten(repmat(uniqueGroups, nReps, 1))';
else
    groups = uniqueGroups;
end

end