function [group, shuffle] = chainGroupTransformations(varargin)
% chainGroupTransformations Perform a series of symmetry transformations on
% a texture group.
%   [group, shuffle] = chainGroupTransformations(trafo1, ..., trafoN, group)
%   applies the transformations in order from left to right (first `trafo1`
%   then `trafo2`, etc.) on the given texture `group`, and composes the
%   returned `shuffle` values.

shuffle = [];
group = varargin{end};
for i = 1:length(varargin)-1
    trafo = varargin{i};
    [group, reshuffle] = trafo(group);
    if isempty(shuffle)
        shuffle = reshuffle;
    else
        shuffle = shuffle(reshuffle);
    end
end

end