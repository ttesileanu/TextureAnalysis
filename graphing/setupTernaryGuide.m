function setupTernaryGuide(groupInfo, varargin)
% setupTernaryGuide Set up the axes with appropriate guides for
% representing a ternary texture plane.
%   setupTernaryGuide(nGroups) prepares the current axes for representing a
%   ternary texture plane with `nGroups` group. `nGroups` can be 1 or 2.
%   This sets `axis` to 'equal' and draws the appropriate guides.
%
%   setupTernaryGuide(groupName) uses the group name to infer the number of
%   groups.
%
%   Options:
%    'triangleOptions'
%       Options to pass to `drawTernaryTriangle`.
%    'mixedOptions'
%       Options to pass to `drawTernaryMixedBackground`.
%   See also: drawTernaryTriangle, drawTernaryMixedBackground.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('triangleOptions', {}, @(c) isempty(c) || (iscell(c) && isvector(c)));
parser.addParameter('mixedOptions', {}, @(c) isempty(c) || (iscell(c) && isvector(c)));
parser.addParameter('labelOptions', {}, @(c) isempty(c) || (iscell(c) && isvector(c)));

% show defaults if requested
if nargin == 1 && strcmp(groupInfo, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse options
parser.parse(varargin{:});
params = parser.Results;

% handle the two kinds of input
if ischar(groupInfo)
    nGroups = 1 + sum(groupInfo == ';');
else
    nGroups = groupInfo;
end

% draw guides
switch nGroups
    case 1
        drawTernaryTriangle(params.triangleOptions{:});
    case 2
        drawTernaryMixedBackground(params.mixedOptions{:});
    otherwise
        error([mfilename ':badngrp'], 'This function only works with single groups and pairs of groups.');
end

axis equal;

end
