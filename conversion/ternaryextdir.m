function v = ternaryextdir(group, direction)
% TERNARYEXTDIR Extend {group, direction} pair(s) to (a) 99-dimensional
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
%
% See also: TERNARYCOMPRESS.

% handle cell array input
if iscell(group) || iscell(direction)
    if ~iscell(group) || ~iscell(direction) || ~isequal(size(group) , size(direction))
        error([mfilename ':badinp'], 'The input should have the same number of groups and directions.');
    end
    v = cellfun(@(g, d) ternaryextdir(g, d), group, direction, 'uniform', false);
    return;
end

% mapping from the 99 probabilities to texture groups
coordGroups = getCoordinateMapping(3, 'expanded');

% all coordinates are centered unless explicitly changed
v = 1/3*ones(1, length(coordGroups));

% update coordinates group by group
sepGroups = strtrim(strsplit(group, ';'));
for k = 1:length(sepGroups)
    crtGroup = sepGroups{k};
    % get rid of a direction specification (in square brackets), if present
    % (that's already included in the direction)
    sqbrIdx = find(crtGroup == '[', 1);
    if ~isempty(sqbrIdx)
        crtGroup = crtGroup(1:sqbrIdx-1);
    end
    crtDirection = direction(3*k-2:3*k);
    crtMask = strcmp(coordGroups, crtGroup);
    v(crtMask) = crtDirection;
end

end