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
%   Options:
%    'covariance': logical
%       Set to true (default) to evaluate a covariance matrix for the
%       patches.
%    'preprocessing':
%       A cell array of preprocessing functions that will be applied to
%       each image, in the same format as for walkImages.
%    Any other key-value arguments are directly passed to walkImages.
%
%   See also: walkImages, analyzeTexture, patchify, getImageTextures.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
parser.KeepUnmatched = true;

parser.addParameter('covariance', true, @(b) islogical(b) && isscalar(b));
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
results = walkImages([params.preprocessing {@walker}], images, unmatched);

% generate covariance, if requested
if params.covariance
    results.covM = safeCov(results.ev);
end

    function crtRes = walker(i, image, crop)
        % walker(i, image, crop) processes the image and returns texture
        % statistics data for it.
        
        % calculate statistics
        try
            crtRes = getImageTextures(image, nLevels);
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
        crtRes.patchLocations = crop;
        crtRes.imageIds = i*ones(size(crtRes.ev, 1), 1);
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