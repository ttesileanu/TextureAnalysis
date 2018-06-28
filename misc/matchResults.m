function [idxs1, idxs2] = matchResults(results1, results2, varargin)
% matchResults Find a mapping between two sets of texture results obtained
% for the same set of images.
%   [idxs1, idxs2] = matchResults(results1, results2) finds matching
%   patches between the two sets of results and returns two vectors of
%   indices identifying the matching patches (i.e., for every `i`, the
%   patch at position `idxs1(i)` in `results1` matches the patch at
%   position `idxs2(i)` in `results`).
%
%   The match uses the image IDs and the patch locations, but does not
%   check whether the image names match. Note that this function is not
%   particularly efficient and can be pretty slow.
%
%   Options:
%    'idFields' {idField1, idField2}
%       Name of fields to read image IDs from in the two results structures.
%    'locationFields' {locationField1, locationField2}
%       Name of fields to read patch location information from in the two
%       results structures.
%
%   See also: generateTextureDistribution.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('idFields', {'imageIds', 'imageIds'}, @(c) ...
    iscell(c) && isvector(c) && length(c) == 2);
parser.addParameter('locationFields', {'patchLocations', 'patchLocations'}, ...
    @(c) iscell(c) && isvector(c) && length(c) == 2);

% defaults
if nargin == 1 && strcmp(results1, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% generate 'hashes' for each results structure
hashes1 = getHashes(results1, params.idFields{1}, params.locationFields{1});
hashes2 = getHashes(results2, params.idFields{2}, params.locationFields{2});

% use the hashes to make a map from patch identity to index
map1 = containers.Map(hashes1, 1:length(hashes1));
map2 = containers.Map(hashes2, 1:length(hashes2));

% find the intersection of the patches, and return the matching indices
commonHashes = intersect(map1.keys, map2.keys);
idxs1 = cell2mat(map1.values(commonHashes));
idxs2 = cell2mat(map2.values(commonHashes));

end

function hashes = getHashes(results, idField, locationField)
% Generate hashes. Note that this assumes that all image sizes and image
% IDs are less than 100,000.

ids = results.(idField);
locs = results.(locationField);

mat = [ids(:) locs];
mats = num2str(mat, '%05d');

hashes = arrayfun(@(i) mats(i, :), 1:size(mats, 1), 'uniform', false);

end