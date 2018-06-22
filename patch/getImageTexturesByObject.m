function res = getImageTexturesByObject(image, nLevels, mask, varargin)
% getImageTexturesByObject Calculate texture statistics for objects in an
% image.
%   res = getImageTexturesByObject(image, nLevels, mask) calculates texture
%   statistics for each object in the scene. Objects are identified by the
%   unique entries in the matrix `mask`. The texture statistics are
%   calculated using `nLevels` levels. If the `mask` is an empty matrix, it
%   is assumed that there is only one object.
%
%   res = getImageTexturesByObject(image, nLevels, mask, 'patchSize', patchSize)
%   calculates texture statistics by splitting the objects into rectangular
%   patches, and calculating statistics for each patch. The options from
%   `getImageTextures` can be used to tune the way in which patches are
%   generated.
%
%   Options:
%     See options for `getImageTextures`. Note that `maxPatchesPerImage`
%     will apply per image *and per object*.
%
%   The output is a structure with the following fields:
%    'objIds': vector
%       The identity of the object analyzed for each patch.
%    'ev': [nPatches, nStats] matrix
%       Matrix containing the statistics for each of the patches.
%    'pxPerPatch': vector
%       The number of pixels that are analyzed for each patch.
%    'patchLocations': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels relative
%       to `image` coordinates.
%    'patchLocationsOrig': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels relative
%       to `mask` coordinates (see `maskCrop` option).
%
%   See also: getImageTextures.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('patchSize', [], @(v) isempty(v) || (isnumeric(v) && ...
    isvector(v) && ismember(length(v), [1 2])));

% defaults
if nargin == 1 && strcmp(image, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% find all the objects in the mask
if ~isempty(mask)
    allObjIds = unique(mask(:));
else
    allObjIds = 1;
end

% split the scene into objects and calculate statistics for each
ev = [];
pxPerPatch = [];
objIds = [];
locs = [];
locsOrig = [];
for i = 1:numel(allObjIds)
    crtObjId = allObjIds(i);
    
    % focus on a given object
    if ~isempty(mask)
        objMask = (mask == crtObjId);
    else
        objMask = [];
    end
    crtRes = getImageTextures(image, nLevels, params.patchSize, 'mask', objMask, unmatched);
    
    ev = [ev ; crtRes.ev]; %#ok<AGROW>
    pxPerPatch = [pxPerPatch ; crtRes.pxPerPatch]; %#ok<AGROW>
    locs = [locs ; crtRes.patchLocations]; %#ok<AGROW>
    locsOrig = [locsOrig ; crtRes.patchLocationsOrig]; %#ok<AGROW>
    objIds = [objIds ; double(crtObjId)*ones(size(crtRes.ev, 1), 1)]; %#ok<AGROW>
end

% fill the output structure
res = struct;
res.objIds = objIds;
res.ev = ev;
res.pxPerPatch = pxPerPatch;

if ~isempty(params.patchSize)
    res.patchLocations = locs;
    res.patchLocationsOrig = locsOrig;
end

end