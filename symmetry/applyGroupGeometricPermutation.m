function [outGroup, shuffle] = applyGroupGeometricPermutation(group, G, permutation)
% applyGroupGeometricPermutation Calculate the effect of a geometric
% permutation on a texture group.
%   [outGroup, permutation] = applyGroupGeometricPermutation(group, G, permutation)
%   calculates how a permutation of the A, B, C, D components within a 2x2
%   grid affects a texture group with `G` gray levels. In more detail, the
%   `permutation` should be a string showing how the values in each 2x2
%   block,
%       A B
%       C D
%   change under the transformation. This is meant to be applied for
%   geometric transformations (e.g., a horizontal flip exchanges A with B
%   and C with D, and would be represented by `permutation = 'BADC'`;
%   clockwise rotation by 90 degress would be `permutation = 'BDAC'`), but
%   can in principle also model transformations that do not map to any
%   global geometric change. This will in general also shuffle the
%   coordinates within the group; the specific shuffling is returned in
%   the matix `shuffle`. Specifically, a direction vector `v` in the
%   initial group gets transformed into direction `v*shuffle` in the
%   `outGroup`. Here `v` is assumed to be a row vector.
%
%   The function also works for mixed groups, separated by a semicolon in
%   the `group` argument.
%
%   applyGroupGeometricPermutation({group1, group2, ...}, ...) applies the
%   transformation to each group, returning cell arrays for `outGroup` and
%   `shuffle`.
%
%   Note that currently this function only works if `G` is a prime number!

% make sure we're using a prime G
if ~isprime(G)
    error([mfilename ':notimp'], 'Currently only prime values of G are supported.');
end

% handle cell array input
if iscell(group)
    [outGroup, shuffle] = cellfun(...
        @(g) applyGroupColorTransformation(g, G, x, y), group, 'uniform', false);
    return;
end

% find out what kind of plane we're looking at
groups = strtrim(strsplit(group, ';'));

% the directions in the transformed groups might get reshuffled
% shuffle0 = zeros(1, G*length(groups));
shuffle0 = zeros(G*length(groups));

% perform the transformation group by group
outGroups0 = cell(size(groups));

% we do so group-by-group
for i = 1:length(groups)
    % parse the group name
    crtGroup = groups{i};
    
    % split out a potential direction specification
    bracketIdx = find(crtGroup == '[');
    if ~isempty(bracketIdx)
        crtDirection = str2double(crtGroup(bracketIdx+1));
        crtGroup = crtGroup(1:bracketIdx-1);
    else
        crtDirection = [];
    end
    
    groupParts = strsplit(crtGroup, '_');
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
    
    % for group orders <= 2, stationarity implies that some groups are equivalent
    % thus normalize the new group name to a representative for each equivalence class
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
    outGroup = cell2mat(newGroupParts);
    
    % figure out how the coordinates within the original group map to
    % coordinates within the target group
    % the mapping is non-trivial when scaling ~= 1
    
    % a combination with index p in the original group will be indexed
    % coordinateMapping(p) in the final group
    coordinateMapping = 1 + mod(scaling*(0:G-1), G);
    
    if ~isempty(crtDirection)
        outDirection = coordinateMapping(crtDirection + 1) - 1;
        outGroup = [outGroup '[' int2str(outDirection) ']']; %#ok<AGROW>
    end
    
    outGroups0{i} = outGroup;
    
    % fill out the shuffle vector
    groupIdxRange = G*(i-1)+1:G*i;
%     shuffle0(groupIdxRange(coordinateMapping)) = groupIdxRange;
    shuffle0(groupIdxRange(coordinateMapping), groupIdxRange) = eye(G);
end

% use canonical ordering of group names
[outGroups, outOrdering] = sortGroups(outGroups0);
shuffle = zeros(size(shuffle0));
for i = 1:length(outGroups)
    j = outOrdering(i);
    shuffle(:, 3*i-2:3*i) = shuffle0(:, 3*j-2:3*j);
end

% reconstruct group name
outGroup = buildGroupName(outGroups{:});

end

function y = modInv(x, G)
% Modular inverse.

[~, a] = gcd(x, G);
y = mod(a, G);

end