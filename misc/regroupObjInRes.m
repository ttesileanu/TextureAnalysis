function newRes = regroupObjInRes(res, labelGroups, varargin)
% regroupObjInRes Regroup object IDs in a texture statistics structure.
%   newRes = regroupObjInRes(res, {labels1, labels2, ...}) returns a new
%   statistics structure in which the object IDs from `labels1` are
%   replaced by 1, the ones in `labels2` are replaced by 2, etc. Patches
%   that don't appear in any of the label groups are removed from the
%   structure.
%
%   regroupObjInRes(res, {labels1, labels2, ...}, [newLabel1 newLabel2 ...])
%   uses `newLabel1`, `newLabel2`, etc. instead of 1, 2, ... for the new
%   object labels.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('newLabels', [], @(v) isvector(v) && isnumeric(v));

% parse
parser.parse(varargin{:});
params = parser.Results;

% set default new labels
if isempty(params.newLabels)
    params.newLabels = 1:length(labelGroups);
end

% update object IDs
newRes = res;
newRes.objIds = nan(size(newRes.objIds));
for i = 1:length(params.newLabels)
    newRes.objIds(ismember(res.objIds, labelGroups{i})) = params.newLabels(i);
end

% remove patches with labels that did not make it into the label groups
newRes = selectRes(newRes, ~isnan(newRes.objIds));

end