function [transformed, shuffle] = applyToThresholds(data, fct, varargin)
% applyToThresholds Apply a function to a structure containing
% psychophysical thresholds.
%   transformed = applyToThresholds(data, fct) applies a transformation
%   function `fct` to each threshold measurement from the structure `data`
%   (as returned by `loadTernaryPP`, for instance), returning the resulting
%   structure with transformed threshold directions. By default, the
%   ordering of the data is unchanged. Use the 'closed' option (below) to
%   restrict the output to those directions that were contained in the
%   initial data. The function needs `data` to have at least `groups` and
%   `directions` fields. The transformation function `fct` needs to take in
%   a group name, and return a transformed group and a shuffling vector.
%   A field called `changed` will be added to the output: it is a logical
%   vector indicating which of the entries were changed by the
%   transformation, and which were left invariant.
%
%   [transformed, shuffle] = applyToThresholds(data, fct) returns a vector
%   `shuffle` showing how the transformation turned one direction into
%   another in the dataset. Specifically, `shuffle(i)` is the index in
%   `data` corresponding to the transformed direction. If the transformed
%   direction is not contained in the dataset, `shuffle(i)` is zero. Note
%   that the 'closed' option (below) effectively only keeps the entries for
%   which `shuffle` is nonzero. Note also that the `changed` field above is
%   simply given by `shuffle(i) ~= i`.
%
%   If the 'closed' option is not used, only the `groups` and `directions`
%   fields are changed in the output. When 'closed' is set to true, all
%   other members of the `transformed` structure that have the same length
%   (if they are vectors) or number of rows (if they are matrices) are
%   subselected in the same way as `groups` and `directions`.
%
%   Options:
%    'closed'
%       When set to true, only those entries that, after the transformation,
%       have directions that were contained in the initial `data` are
%       returned.

% setup the input options
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('closed', false, @(b) isscalar(b) && islogical(b));

% show defaults
if nargin == 1 && strcmp(data, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% apply transformation
transformed = data;
shuffle = zeros(1, length(data.groups));
for i = 1:length(data.groups)
    [trafoGroup, crtShuffle] = fct(data.groups{i});
    trafoDir = data.directions{i}(crtShuffle);
    
    % store the transformed data
    transformed.groups{i} = trafoGroup;
    transformed.directions{i} = trafoDir;
    
    % is the transformed direction in the data?
    groupMatchMask = strcmp(data.groups, trafoGroup);
    groupMatchIdxs = find(groupMatchMask);
    if ~isempty(groupMatchIdxs)
        % at least the transformed group is...
        directionMatchSubmask = cellfun(@(d) max(abs(d - trafoDir)) < 1e-6, ...
            data.directions(groupMatchMask));
        nMatches = sum(directionMatchSubmask);
        % the direction, too
        if nMatches > 0
            if nMatches > 1
                % but it appears in several locations!
                warning([mfilename ':multimatch'], 'Data at index %d (group %s) matches several datapoints after transformation.', ...
                    i, data.groups{i});
            end
            submaskIdx = find(directionMatchSubmask, 1);
            shuffle(i) = groupMatchIdxs(submaskIdx);
        end
    end
end

% find which directions changed
transformed.changed = flatten(shuffle ~= 1:length(shuffle));

% if asked to, keep only directions that existed in original data
if params.closed
    transformed = selectMeasurements(transformed, shuffle > 0);
end

end