function group = normalizeGroup(group)
% normalizeGroup Normalize texture group name.
%   group = normalizeGroup(group) normalizes the group name by, e.g.,
%   padding it with spaces so that it sorts in a natural way and can be
%   used for comparisons.
%
%   If `group` is a cell array, the function normalizes each element. The
%   function also works for multi-group names (where the groups are
%   separated by semicolons).

% handle empty input
if isempty(group)
    return;
end

% handle cell array version
if iscell(group)
    group = cellfun(@normalizeGroup, group, 'uniform', false);
    return;
end

% handle semicolon-separated multi-group names
multiGroup = strtrim(strsplit(group, ';'));
if length(multiGroup) > 1
    multiGroup = normalizeGroup(multiGroup);
    parts = cell(1, 2*length(multiGroup)-1);
    parts(1:2:end) = multiGroup;
    parts(2:2:end) = {';'};
    group = strcat(parts{:});
    return;
end

% remove white space
group = group(~isspace(group));

underStart = find(group == '_', 1);
letters0 = group(1:underStart-1);
letters = pad(letters0, 4, 'left');

dirStart = find(group == '[', 1);
if ~isempty(dirStart)
    dirEnd = find(group == ']', 1);
    direction = group(dirStart:dirEnd);
else
    dirStart = length(group) + 1;
    direction = '[0]';
end

under0 = group(underStart:dirStart-1);
under = [repmat('_0', 1, 4-length(letters0)) under0];

group = [letters under direction];

end