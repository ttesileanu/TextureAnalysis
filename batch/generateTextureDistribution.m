function results = generateTextureDistribution(images, nLevels, patchSize, varargin)
% generateTextureDistribution Calculate texture statistics for a set of
% images.
%   results = generateTextureDistribution(images, nLevels, patchSize)
%   calculates texture statistics with `nLevels` gray levels and the given
%   `patchSize` for the images in `images` (which can be in any of the
%   formats accepted by walkImages). Note that the texture analysis routine
%   assumes that the preprocessing routines (see 'preprocessing' option
%   below) convert the image to a format in which all pixels take `nLevels`
%   discrete values in the interval [0, 1] (except in the case in which
%   `nLevels` is infinite; see analyzeTexture).
%
%   results = generateTextureDistribution(images, 'masks', masks) uses the
%   images given in the cell array `masks` to identify portions of the
%   image on which to focus the analysis (see XXX for details).
%   Empty entries in the `masks` cell array can be used to skip some of the
%   images in the image set.
%
%   Options:
%    'stride': scalar or pair of numbers
%       If set, use a sliding patch instead of being restricted to
%       overlapping patches. See `ImagePatchifier`.
%    'masks': cell array of matrices
%       The images in this cell array are used to identify portions of the
%       image on which to focus the analysis (see XXX for details).
%       Empty entries in the `masks` cell array can be used to skip some of
%       the images in the image set.
%    'maxPatchesPerImage': integer
%       If given, this ensures that the number of patches per image doesn't
%       exceed a given value. If an image generates more patches than
%       `maxPatchesPerImage`, random sampling without replacement is used
%       to select `maxPatchesPerImage` of them to keep.
%       (default: no maximum)
%    'minPatchUsed': double
%       Minimum fraction of patch that needs to be available (e.g.,
%       contained within a region identified by the masks) to be analyzed.
%       See XXX.
%    'covariances': logical
%       Set to true (default) to evaluate covariance matrices for the
%       patches. If no masks are provided, a single covarince matrix is
%       calculated, averaging over all patches. If there are masks, then a
%       covariance matrix is calculated for every object ID found in the
%       masks.
%    'preprocessing':
%       A cell array of preprocessing functions that will be applied to
%       each image, in the same format as for walkImages.
%    Any other key-value arguments are directly passed to walkImages.
%
%   See also: walkImages, analyzeTexture, ImagePatchifier.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('masks', {}, @(c) iscell(c));
parser.addParameter('stride', [], @(v) isempty(v) || ...
    (isnumeric(v) && isvector(v) && ismember(length(v), [1, 2]) && all(v > 0)));
parser.addParameter('maxPatchesPerImage', inf, @(n) isnumeric(n) && isscalar(n) && n >= 0);
parser.addParameter('minPatchUsed', 0, @(n) isnumeric(n) && isscalar(n) && n >= 0 && n <= 1);
parser.addParameter('covariances', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('preprocessing', {}, @(c) iscell(c) && (isempty(c) || isvector(c)));

% defaults
if nargin == 1 && strcmp(images, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;
unmatched = parser.Unmatched;

% set up the arguments to pass to getImageTexturesByObject
analysisArgs = {'patchSize', patchSize, 'stride', params.stride, ...
    'maxPatchesPerImage', params.maxPatchesPerImage, 'minPatchUsed', params.minPatchUsed};

% process the images
masks = params.masks;
results = walkImages([params.preprocessing {@walker}], images, unmatched);

% generate covariances, if requested
if params.covariances
    results.covM = safeCov(results.ev);
    if ~isempty(params.masks)
        % this assumes that object IDs are consistent across images!
        allObjs = unique(results.objIds);
        results.covPerObj = cellfun(@(crtObj) ...
            safeCov(results.ev(results.objIds == crtObj)), 1:length(allObjs));
    end
end

% keep track of the options used to generate the texture distribution
results.options.patchSize = patchSize;
results.options.stride = params.stride;
results.options.maxPatchesPerImage = params.maxPatchesPerImage;
results.options.minPatchUsed = params.minPatchUsed;

    function crtRes = walker(i, image, crop)
        % walker(i, image, crop) processes the image and returns texture
        % statistics data for it.
        
        if isempty(masks)
            % there is no mask
            mask = [];
        else
            % there are masks
            mask = masks{i};
            % if this particular mask is empty, skip image
            if isempty(mask)
                crtRes = [];
                return;
            end
        end
        
        % calculate statistics
        try
            crtRes = getImageTexturesByObject(image, nLevels, mask, analysisArgs{:}, ...
                'maskCrop', crop);
        catch errMsg
            newErrMsg.message = ['Image #' int2str(i) ': ' errMsg.message];
            newErrMsg.identifier = errMsg.identifier;
            newErrMsg.stack = errMsg.stack;
            error(newErrMsg);
        end
        if isempty(crtRes.ev)
            % skip images for which we got no texture data
            crtRes = [];
            return;
        end
        % keep track of the image from which the patches originated
        crtRes.imgIds = i*ones(size(crtRes.ev, 1), 1);
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