function [ev, shuffle] = applyStatsGeometricPermutation(ev, G, permutation)
% applyStatsGeometricPermutation Calculate the effect of a geometric
% permutation on texture statistics.
%   evP = applyStatsGeometricPermutation(ev, G, perumtation) calculates the
%   way in which texture statistics with `G` gray levels contained in the
%   matrix `ev` would change if the patches from which these statistics were
%   calculated were modified by a permutation of the A, B, C, D components
%   within a 2x2 grid. In more detail, the `permutation` should be a string
%   showing how the values in each 2x2 block,
%       A B
%       C D
%   change under the transformation. This is meant to be applied for
%   geometric transformations (e.g., a horizontal flip exchanges A with B
%   and C with D, and would be represented by `permutation = 'BADC'`;
%   clockwise rotation by 90 degress would be `permutation = 'BDAC'`), but
%   can in principle also model transformations that do not map to any
%   global geometric change.
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
% hard-coding the mapping for binary and ternary textures to reduce
% dependence on Jonathan's framework
if G == 2
    coordGroups = {'A_1'    'AB_1_1'    'AC_1_1'    'BC_1_1'    'AD_1_1' ...
        'ABC_1_1_1'    'ABD_1_1_1'    'ACD_1_1_1'    'BCD_1_1_1'    'ABCD_1_1_1_1'};
elseif G == 3
    coordGroups = {...
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
    coordGroups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
end

% the transformed stats will be related to the original ones by a shuffling
% of columns; we aim to find the shuffling vector
shuffle = zeros(1, size(ev, 2));

% we do so group-by-group
for i = 1:length(coordGroups)
    % parse the group name
    group = coordGroups{i};
    groupParts = strsplit(group, '_');
    groupLetters = groupParts{1};
    groupCoefficients = str2double(groupParts(2:end));
    
    % permute the letters
    % XXX this is a bit hacky
    groupLettersAsNum = groupLetters - 'A' + 1;
    newGroupLetters = permutation(groupLettersAsNum);
    
    % now normalize the group name again, keeping consistency with the
    % coefficients
    [newGroupLetters, reordering] = sort(newGroupLetters);
    newGroupCoefficients = groupCoefficients(reordering);
    
    % for group orders up to 2, stationarity implies that some groups are
    % equivalent
    % we thus normalize the new group name to a representative for each
    % equivalence class
    if isscalar(newGroupLetters)
        % the only first-order statistic corresponds to 'A'
        newGroupLetters = 'A';
    elseif length(newGroupLetters) == 2
        switch newGroupLetters
            case 'CD'
                newGroupLetters = 'AB';
            case 'BD'
                newGroupLetters = 'AC';
        end
        % other cases are the representatives of their own classes
    end
    
    % now make sure the first coefficient is 1
    scaling = modInv(newGroupCoefficients(1), G);
    newGroupCoefficients = mod(scaling*newGroupCoefficients, G);
    
    % generate the name of the group to which this group maps
    newGroupParts = [{newGroupLetters} ...
        flatten([repmat({'_'}, 1, length(newGroupCoefficients)) ; ...
                 arrayfun(@int2str, newGroupCoefficients, 'uniform', false)])'];
    newGroup = cell2mat(newGroupParts);
    
    % figure out how the coordinates within the original group map to
    % coordinates within the target group
    % the mapping is non-trivial when scaling ~= 1
    
    % a combination with index p in the original group will be indexed
    % coordinateMapping(p) in the final group
    coordinateMapping = 1 + mod(scaling*(0:G-1), G);
    
    % fill out the shuffle vector
    newGroupIdx = find(strcmp(coordGroups, newGroup));
    groupIdxRange = G*(i-1)+1:G*i;
    newGroupIdxRange = G*(newGroupIdx-1)+1:G*newGroupIdx;
    shuffle(newGroupIdxRange(coordinateMapping)) = groupIdxRange;
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