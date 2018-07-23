function group = buildGroupName(varargin)
% buildGroupName Assemble a mixed group name from a sequence of groups.
%   group = buildGroupName(group1, group2, ...) creates a single group name
%   by concatenating the given groups, separated by semicolons.

group = cell2mat(flatten([varargin ; [repmat({';'}, 1, length(varargin)-1) {''}]])');

end