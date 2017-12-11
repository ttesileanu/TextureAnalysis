function res = analyzeImageSet(imageNames, path, varargin)
% analyzeImageSet Calculate texture statistics for a set of images.
%   res = analyzeImageSet(imageNames, path) calculates continuous texture
%   statistics for the images located in folder `path` with names given by
%   the `imageNames` cell array. See the options below (and in
%   `walkImageSet`) for preprocessing options.
%
%   res = analyzeImageSet(imageNames, path, blockAF) performs block
%   averaging by the given factor before running the analysis.
%
%   res = analyzeImageSet(..., 'nLevels', n) performs an `n`-graylevel
%   analysis instead of the continuous texture analysis. This is implied by
%   default when the `quantize` option is used. Conversely, if `nLevels` is
%   provided, the quantization is implied.
%
%   res = analyzeImageSet(imageNames, path, masks, ...) uses the images
%   given in the cell array `masks` to identify portions of the image on
%   which to focus the analysis (see `analyzeObjects` for details). Empty
%   entries in the `masks` cell array can be used to skip some of the
%   images in the image set.
%
%   res = analyzeImageSet(imageCount, imageGenerator, ...) uses a function
%   handle (`imageGenerator`) to generate the `imageCount` images that are
%   to be processed. `imageGenerator` is called as `imageGenerator(i)`,
%   where `i` is the index of the image that's being requested, and should
%   return an image that can be handled by `preprocessImage`.
%
%   Options:
%    'averageType':
%    'doLog':
%    'threshold':
%    'filter':
%    'filterType':
%    'equalize':
%    'equalizeType': 
%    'patchSize':
%    'nonlinearity':
%    'quantize':
%    'quantizeType':
%       These control the preprocessing of the images. See `walkImageSet`
%       for a description. Note that `patchSize` is by default set to
%       the `size(filter)`, if the `filter` is not empty, or to
%       `analysisPatchSize` (if provided).
%    'analysisPatchSize': integer, or pair of integers
%       When provided, the analysis is performed by splitting each image
%       (or each region within an image as identified by the masks) into
%       patches of the given size. This can be a single number to use
%       square patches, or a pair of numbers [sizeY, sizeX] for rectangular
%       patches. If 'analysisPatchSize' is not provided, it is inferred
%       from 'filter' or 'patchSize', if they are provided. To do the
%       texture analysis on the entire image regardless of the size of the
%       whitening filter or of the patch size used for
%       equalization/quantization, set 'analysisPatchSize' to an empty
%       matrix.
%    'overlapping': logical, numeric, or pair of numbers
%       Set to true, or to a step size, or to a pair of step sizes, to use
%       a pixel-by-pixel sliding patch to evaluate statistics. The patches
%       are centered on a pixel that's sliding with the given step. This
%       can generate many more patches but the patches are no longer
%       independent of each other. See `analyzePatches` and
%       `analyzeObjects`. Note that this option is ignored if no
%       `patchSize` is given.
%    'maxPatchesPerImage': integer
%       If given, this ensures that the number of patches per image doesn't
%       exceed a given value. If an image generates more patches than
%       `maxPatchesPerImage`, random sampling without replacement is used
%       to select `maxPatchesPerImage` of them to keep.
%       (default: no maximum)
%    'minPatchUsed': double
%       Minimum fraction of patch that needs to be available (e.g.,
%       contained within a region identified by the masks) to be analyzed.
%       See `analyzePatches`.
%    'nLevels': int
%       Number of graylevels to use for the texture analysis. This
%       overrides the 'quantize' option.
%       (default: the value from 'quantize' if provided, otherwise use
%                 continuous analysis)
%    'covariances': logical
%       Set to true (default) to evaluate covariance matrices for the
%       patches. If no masks are provided, a single covarince matrix is
%       calculated, averaging over all patches. If there are masks, then a
%       covariance matrix is calculated for every object ID found in the
%       masks.
%    'imageCopies': logical
%       If true, copies of the original, log-transformed, block-averaged,
%       filtered, and final, quantized, images are included in the output
%       structure. The (downsampled) masks are also included. This is true
%       by default for at most 10 images, false otherwise.
%    'progressEvery':
%    'progressStart':
%       See `walkImageSet`.
%
%   See also: walkImageSet, preprocessImage, analyzeObjects, analyzePatches.

% handle masks
if iscell(varargin{1})
    masks = varargin{1};
    varargin = varargin(2:end);
else
    masks = {};
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

checkStr = @(s) isempty(s) || (ischar(s) && isvector(s));
checkBool = @(b) isempty(b) || (islogical(b) && isscalar(b));
checkPatchSize = @(v) isempty(v) || (isnumeric(v) && isvector(v) && ...
    (numel(v) == 1 || numel(v) == 2) && all(v >= 1));
checkNumber = @(x) isempty(x) || (isscalar(x) && isreal(x) && isnumeric(x));

parser.addOptional('blockAF', 1, checkNumber);

parser.addParameter('nLevels', [], checkNumber);

parser.addParameter('averageType', [], checkStr);
parser.addParameter('doLog', [], checkBool);
parser.addParameter('threshold', [], checkNumber);
parser.addParameter('filter', [], @(m) isempty(m) || (ismatrix(m) && isreal(m) && isnumeric(m)));
parser.addParameter('filterType', [], checkStr);
parser.addParameter('equalize', [], checkBool);
parser.addParameter('equalizeType', [], checkStr);
parser.addParameter('patchSize', [], checkPatchSize);
parser.addParameter('nonlinearity', [], @(v) isempty(v) || (isvector(v) && isnumeric(v) && isreal(v)));
parser.addParameter('quantize', [], checkNumber);
parser.addParameter('quantizeType', [], checkStr);
parser.addParameter('analysisPatchSize', 'default', checkPatchSize);

parser.addParameter('overlapping', []);
parser.addParameter('maxPatchesPerImage', inf, checkNumber);
parser.addParameter('minPatchUsed', [], checkNumber);
parser.addParameter('covariances', true, checkBool);

parser.addParameter('imageCopies', [], checkBool);
parser.addParameter('progressEvery', [], checkNumber);
parser.addParameter('progressStart', [], checkNumber);

% parse
parser.parse(varargin{:});
params = parser.Results;

% default patchSize is the same as size of filter, or as analysisPatchSize
if isempty(params.patchSize)
    if isempty(params.filter)
        if ~strcmp(params.analysisPatchSize, 'default')
            params.patchSize = params.analysisPatchSize;
        else
            params.patchSize = [];
        end
    else
        params.patchSize = size(params.filter);
    end
end

% default analysisPatchSize is the same as patchSize or filter
if strcmp(params.analysisPatchSize, 'default')
    if ~isempty(params.filter)
        params.analysisPatchSize = size(params.filter);
    elseif ~isempty(params.patchSize)
        params.analysisPatchSize = params.patchSize;
    end
end

% set defaults for nLevels and quantize
if ~isempty(params.nLevels) && ~isempty(params.quantize)
    if isfinite(params.nLevels) && ~isequal(params.nLevels, params.quantize)
        error([mfilename ':badnlev'], 'nLevels and quantize must match when they are both provided, unless nLevels is infinity.');
    end    
elseif ~isempty(params.nLevels)
    params.quantize = params.nLevels;
elseif ~isempty(params.quantize)
    params.nLevels = params.quantize;
end
if isempty(params.nLevels)
    params.nLevels = inf;
end

% set the default for image copies
if isempty(params.imageCopies)
    if iscell(imageNames)
        params.imageCopies = numel(imageNames) <= 10;
    else
        params.imageCopies = false;
    end
end

% generate argument lists for walkImageSet and analyzeObjects
optionalArgs = structToCell(params, ...
    {'averageType', 'doLog', 'equalize', 'equalizeType', 'filter', 'filterType', ...
    'nonlinearity', 'quantize', 'quantizeType', 'patchSize', 'threshold', ...
    'progressEvery', 'progressStart'});
if ~isempty(params.patchSize)
    analysisArgs = {params.analysisPatchSize};
else
    analysisArgs = {};
end
analysisArgs = [analysisArgs structToCell(params, {'minPatchUsed', 'overlapping', 'maxPatchesPerImage'})];

res = walkImageSet(@walker, imageNames, path, params.blockAF, optionalArgs{:});

if params.covariances
    res.covM = safeCov(res.ev);
    if ~isempty(masks)
        % this assumes that object IDs are consistent across images!
        allObjs = unique(res.objIds);
        res.covPerObj = cell(1, length(allObjs));
        for i = 1:length(allObjs)
            crtObj = allObjs(i);
            crtEv = res.ev(res.objIds == crtObj, :);
            res.covPerObj{i} = safeCov(crtEv);
        end
    end
end

res.options.analysisPatchSize = params.analysisPatchSize;
res.options.overlapping = params.overlapping;
res.options.maxPatchesPerImage = params.maxPatchesPerImage;
res.options.minPatchUsed = params.minPatchUsed;

    function crtRes = walker(i, image, details)
        % walker(i, image, details) processes the image and returns texture
        % statistics data for it.
        
        if isempty(image)
            % if preprocessing reduced the image to nothing, skip it
            warning([mfilename ':smallimg'], ['Image #' int2str(i) ' was trimmed to zero size by preprocessing -- skipping.']);
            crtRes = [];
            return;
        end
        
        if isempty(masks)
            % need to create a mask as large as the unscaled & untrimmed
            % image
            mask = ones(details.crops.final(3), details.crops.final(4));
        else
            mask = masks{i};
            % if there is no mask, skip image
            if isempty(mask)
                crtRes = [];
                return;
            end
        end
        
        % calculate statistics
        try
            crtRes = analyzeObjects(image, params.nLevels, mask, analysisArgs{:}, ...
                'maskCrop', details.crops.final);
        catch errMsg
            newErrMsg.message = ['Image #' int2str(i) ': ' errMsg.message];
            newErrMsg.identifier = errMsg.identifier;
            newErrMsg.stack = errMsg.stack;
            error(newErrMsg);
        end
        if isempty(crtRes.ev)
            crtRes = [];
            return;
        end
        crtRes.imgIds = i*ones(size(crtRes.ev, 1), 1);
        
        crtRes = rmfield(crtRes, intersect(...
            {'nLevels', 'patchSize', 'overlapping', 'minPatchUsed', 'maxPatchesPerImage'}, ...
            fieldnames(crtRes)));
        
        % keep track of the source images if asked to do so
        if params.imageCopies
            crtRes.imageCopies.original = details.origImage;
            crtRes.imageCopies.log = details.logImage;
            crtRes.imageCopies.blockAveraged = details.averagedImage;
            crtRes.imageCopies.filtered = details.filteredImage;
            crtRes.imageCopies.final = image;
            crtRes.imageCopies.mask = mask;
            crtRes.imageCopies.crops = details.crops;
        end
    end

end

function c = safeCov(m)
% safeCov(m) returns the covariance matrix for m, or a matrix of NaN is m
% has only one row. An empty matrix is returned if m is empty.

if isempty(m)
    c = [];
elseif size(m, 1) == 1
    c = nan(size(m, 2));
else
    c = cov(m);
end

end