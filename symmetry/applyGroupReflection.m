function [group, shuffle] = applyGroupReflection(group, G, varargin)
% applyGroupReflection Perform a reflection around the origin of a texture
% group.
%   [group, shuffle] = applyGroupReflection(group, G) calculates the effect
%   of a reflection around the origin (probabilities equal to `1/G`) in the
%   given texture group. The `group` is returned for compatibility with
%   other symmetry functions, but is always unchanged. The transformation
%   is implemented in the `shuffle` matrix, which is defined such as that
%   an input vector `v` gets transformed into `v*shuffle` (with `v` being a
%   row vector).
%
%   If `group` is a cell array, the return values will be cell arrays.
%
%   Options:
%    'restrict': callable
%       If given, the reflection is only performed in specific groups for
%       which the function returns `true`. Other groups are left unchanged.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('restrict', @(g) true);

% show defaults if requested
if nargin == 1 && strcmp(group, 'defaults')
    parser.parse;
    disp(parser.Results);
    group = [];
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% handle cell array input
if iscell(group)
    [group, shuffle] = cellfun(...
        @(g) applyGroupReflection(g, G, varargin{:}), group, 'uniform', false);
    return;
end

% find out what kind of plane we're looking at
groups = strtrim(strsplit(group, ';'));

% the directions will get transformed
shuffle = zeros(G*length(groups));

% perform the transformation group by group
for i = 1:length(groups)
    % parse the group name
    crtGroup = groups{i};
    
    % split out a potential direction specification
    bracketIdx = find(crtGroup == '[');
    if ~isempty(bracketIdx)
        crtGroup = crtGroup(1:bracketIdx-1);
    end

    % fill out the shuffle matrix
    groupIdxRange = G*(i-1)+1:G*i;
    
    % is this group to be transformed?
    if ~params.restrict(crtGroup)
        % no
        shuffle(groupIdxRange, groupIdxRange) = eye(G);
    else
        % yes
        shuffle(groupIdxRange, groupIdxRange) = (2/G)*ones(G) - eye(G);
    end
end

end