function results = generateTextureDistribution(images, nLevels, varargin)
% generateTextureDistribution Calculate texture statistics for a set of
% images.
%   results = generateTextureDistribution(images, nLevels) calculates
%   texture statistics with `nLevels` gray levels for the images in
%   `images` (which can be in any of the formats accepted by walkImages).
%   Note that the texture analysis routine assumes that the preprocessing
%   routines (see 'preprocessing' option below) split the images in
%   appropriately-sized patches (if necessary), and convert them to a
%   format in which all pixels take `nLevels` discrete values in the
%   interval [0, 1] (except in the case in which `nLevels` is infinite; see
%   analyzeTexture).
%
%   results = generateTextureDistribution(..., 'masks', masks) uses the
%   images given in the cell array `masks` to identify portions of each
%   input image on which to focus the analysis (see getImageTexturesByObject
%   for details). Empty entries in the `masks` cell array can be used to
%   skip some of the images in the image set.
%
%   Options:
%    'masks': cell array of matrices
%       The images in this cell array are used to identify portions of the
%       input images on which to focus the analysis (see
%       getImageTexturesByObject for details). Images corresponding to
%       empty entries in the `masks` cell array are skipped.
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
%   See also: walkImages, analyzeTexture, patchify, getImageTextures, getImageTexturesByObject.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('masks', {}, @(c) iscell(c));
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
            crtRes = getImageTexturesByObject(image, nLevels, mask, 'maskCrop', crop);
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