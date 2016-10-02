function res = analyzeObjects(image, nLevels, mask, varargin)
% analyzeObjects Calculate texture statistics for objects in an image.
%   res = analyzeObjects(image, nLevels, mask) calculates texture
%   statistics for each object in the scene. Objects are identified by the
%   unique entries in the matrix `mask`. The texture statistics are
%   calculated using `nLevels` levels.
%
%   res = analyzeObjects(image, nLevels, mask, patchSize) calculates
%   texture statistics by splitting the objects into rectangular patches,
%   and calculating statistics for each patch. The options from
%   `analyzePatches` can be used to tune the way in which patches are
%   generated.
%
%   Options:
%     See options for `analyzePatches`.
%
%   The output is a structure with the following fields:
%    'objIds': vector
%       The identity of the object analyzed for each patch.
%    'ev': [nPatches, nStats] matrix
%       Matrix containing the statistics for each of the patches.
%    'pxPerPatch': vector
%       The number of pixels that are analyzed for each patch.
%    'nLevels':
%       This is just copied from the input argument.
%    'patchLocations': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels relative
%       to `image` coordinates.
%    'patchLocationsOrig': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels relative
%       to `mask` coordinates (see `maskCrop` option).
%    'patchSize':
%    'overlapping':
%    'minPatchUsed':
%       These are just copies of the input arguments and options.
%
%   See also: analyzePatches.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('patchSize', []);

parser.addParameter('minPatchUsed', []);
parser.addParameter('overlapping', []);
parser.addParameter('maskCrop', []);

% parse
parser.parse(varargin{:});
params = parser.Results;

% set up the named arguments for analyzePatches
anPNamedArgs = {};
if ~isempty(params.minPatchUsed)
    anPNamedArgs = [anPNamedArgs {'minPatchUsed' params.minPatchUsed}];
end
if ~isempty(params.overlapping)
    anPNamedArgs = [anPNamedArgs {'overlapping' params.overlapping}];
end
if ~isempty(params.maskCrop)
    anPNamedArgs = [anPNamedArgs {'maskCrop' params.maskCrop}];
end

allObjIds = unique(mask(:));

% split the scene into objects and calculate statistics for each
ev = [];
pxPerPatch = [];
objIds = [];
locs = [];
locsOrig = [];
for i = 1:numel(allObjIds)
    crtObjId = allObjIds(i);
    
    % focus on a given object
    objMask = (mask == crtObjId);
    
    crtRes = analyzePatches(image, nLevels, params.patchSize, objMask, anPNamedArgs{:});
    
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
res.nLevels = nLevels;

if ~isempty(params.patchSize)
    res.patchLocations = locs;
    res.patchLocationsOrig = locsOrig;
    res.patchSize = crtRes.patchSize;
    res.overlapping = crtRes.overlapping;
    res.minPatchUsed = crtRes.minPatchUsed;
end

end