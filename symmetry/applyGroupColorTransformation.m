function [outGroup, shuffle] = applyGroupColorTransformation(group, G, x, y)
% applyGroupColorTransformation Calculate the effect on a texture group of
% applying an affine transformation to gray levels.
%   [outGroup, shuffle] = applyGroupColorTransformation(group, G, x, y)
%   calculates how an affine transformation on the gray levels affects a
%   given texture group. The number of gray levels is assumed to be `G`,
%   and the affine transformation turns level `i` into `mod(x*i + y, G)`.
%   The transformation can turn one texture group into another; the
%   resulting texture group is returned in `outGroup`. The transformation
%   might also result in a shuffling of coordinates within the group, which
%   is returned in the matrix `shuffle`. Specifically, a direction vector
%   `v` in the initial group gets transformed into direction `v*shuffle` in
%   the `outGroup`. Here `v` is assumed to be a row vector.
%
%   The function also works for mixed groups, separated by a semicolon in
%   the `group` argument.
%
%   applyGroupColorTransformation({group1, group2, ...}, ...) applies the
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
% shuffle = zeros(1, G*length(groups));
shuffle = zeros(G*length(groups));

% perform the transformation group by group
outGroups = cell(size(groups));
xInv = modInv(x, G);
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
        
    % the group name stays the same, but the coordinates get shuffled
    % the way in which they do depends on the group
    groupCoefficients = str2double(groupParts(2:end));
    
    % a combination with index p in the original group will be at index
    % coordinateMapping(p) in the final group
    shift = mod(y*sum(groupCoefficients), G);
    coordinateMapping = 1 + mod(xInv*((0:G-1) - shift), G);
    
    % fill out the shuffle vector
    groupIdxRange = G*(i-1)+1:G*i;
%     shuffle(groupIdxRange(coordinateMapping)) = groupIdxRange;
    shuffle(groupIdxRange(coordinateMapping), groupIdxRange) = eye(G);
    
    if ~isempty(crtDirection)
        outDirection = coordinateMapping(crtDirection + 1) - 1;
        crtGroup = [crtGroup '[' int2str(outDirection) ']']; %#ok<AGROW>
    end
    outGroups{i} = crtGroup;
end

% reconstruct group name
outGroup = buildGroupName(outGroups{:});

end

function y = modInv(x, G)
% Modular inverse.

[~, a] = gcd(x, G);
y = mod(a, G);

end