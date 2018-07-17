function [ev, shuffle] = applyStatsColorTransformation(ev, G, x, y)
% applyStatsColorTransformation Calculate the effect on texture statistics
% of applying an affine transformation to gray levels.
%   evP = applyStatsColorTransformation(ev, G, x, y) calculates the way in
%   which texture statistics with `G` gray levels contained in the matrix
%   `ev` would change if the patches from which these statistics were
%   calculated were transformed by an affine transformation of the gray
%   levels, in which level `i` gets transformed into `x*i + y (mod 3)`.
%
%   The result is such that `evP(i, j) == ev(i, shuffle(j))`.
%
%   Note that currently this function only works if `G` is a prime number!
%
%   `ev` can be given either in the independent-component format in which
%   each texture group has `G-1` components, or in the full probability
%   format in which there are `G` probabilities for each texture group.
%   Note that the binary case, `G == 2`, is treated differently, in line
%   with `analyzeTexture`.
%
% See also: analyzeTexture, expandTextureStats.

% make sure we're using a prime G
if ~isprime(G)
    error([mfilename ':notimp'], 'Currently only prime values of G are supported.');
end

% check whether we need to extend stats
nIndependent = G*(G-1)*(G^2 + G - 1);
% keep track of the type in which the initial data was presented, so we can
% revert to it at the end
initialType = 'full';
if ismatrix(ev) && size(ev, 2) == nIndependent
    ev = expandTextureStats(ev, G, '3d');
    initialType = 'independent';
elseif ~ismatrix(ev)
    % transform to a 2d array
    ev = reshape(ev, size(ev, 1), []);
    initialType = '3d';
end

% we need the mapping from groups of G column indices in `ev` vector to
% texture groups
coordGroups = getCoordinateMapping(G);

% the transformed stats will be related to the original ones by a shuffling
% of columns; we aim to find the shuffling vector
shuffle = zeros(1, size(ev, 2));

% we do so group-by-group
xInv = modInv(x, G);
for i = 1:length(coordGroups)
    % parse the group name
    group = coordGroups{i};
    groupParts = strsplit(group, '_');
        
    % the group name stays the same, but the coordinates get shuffled
    % the way in which they do depends on the group
    groupCoefficients = str2double(groupParts(2:end));
    
    % a combination with index p in the original group will be at index
    % coordinateMapping(p) in the final group
    shift = mod(y*sum(groupCoefficients), G);
    coordinateMapping = 1 + mod(xInv*(0:G-1) + shift, G);
    
    % fill out the shuffle vector
    groupIdxRange = G*(i-1)+1:G*i;
    shuffle(groupIdxRange(coordinateMapping)) = groupIdxRange;
end

% shuffle columns of ev
ev = ev(:, shuffle);

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

function y = modInv(x, G)
% Modular inverse.

[~, a] = gcd(x, G);
y = mod(a, G);

end