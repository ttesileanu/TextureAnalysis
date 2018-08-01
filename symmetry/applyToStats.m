function [ev, shuffle] = applyToStats(ev, G, fct)
% applyToStats Find the effect of a symmetry transformation on a matrix of
% texture statistics.
%   evOut = applyToStats(ev, G, fct) calcuales the effect of a symmetry
%   transformation `fct` on texture statistics `ev` with `G` graylevel. The
%   transformation function `fct` takes in a group name and returns a
%   transformed group.
%
%   [evOut, shuffle] = applyToStats(ev, G, fct) returns a matrix `shuffle`
%   showing how to transform the input stats into the output ones,
%       evOut = ev*shuffle .
%
%   `ev` can be given either in the independent-component format in which
%   each texture group has `G-1` components, or in the full probability
%   format in which there are `G` probabilities for each texture group.
%   Note that the binary case, `G == 2`, is treated differently, in line
%   with `analyzeTexture`.
%
%   See also: analyzeTexture, expandTextureStats.

% check whether we need to extend stats
nIndependent = G*(G-1)*(G^2 + G - 1);
% keep track of the type in which the initial data was presented, so we can
% revert to it at the end
initialType = 'full';
if ismatrix(ev) && size(ev, 2) == nIndependent
    ev = expandTextureStats(ev, G);
    initialType = 'independent';
elseif ~ismatrix(ev)
    % transform to a 2d array
    ev = reshape(ev, size(ev, 1), []);
    initialType = '3d';
end

% we need the mapping from groups of G column indices in `ev` vector to
% texture groups
coordGroups = getCoordinateMapping(G);

% find the shuffle matrix
% shuffle = zeros(1, size(ev, 2));
shuffle = zeros(size(ev, 2));

% work group-by-group
for i = 1:length(coordGroups)
    [newGroup, crtShuffle] = fct(coordGroups{i});
    
    % fill out the shuffle vector
    newGroupIdx = find(strcmp(coordGroups, newGroup));
    groupIdxRange = G*(i-1)+1:G*i;
    newGroupIdxRange = G*(newGroupIdx-1)+1:G*newGroupIdx;
%     shuffle(newGroupIdxRange) = groupIdxRange(crtShuffle);
    shuffle(groupIdxRange, newGroupIdxRange) = crtShuffle;
end

% transform ev
% ev = ev(:, shuffle);
ev = ev*shuffle;

% convert back to the initial format for the ev data
switch initialType
    case 'independent'
        % easiest way to do this is to go through the 3d format
        ev = reshape(ev, size(ev, 1), G, []);
        ev = reshape(ev(:, 1:end-1, :), size(ev, 1), []);
    case '3d'
        ev = reshape(ev, size(ev, 1), G, []);
    case 'full'
        % nothing to do
    otherwise
        error([mfilename ':bug'], 'Invalid initialType. This should never hapen.');
end

end