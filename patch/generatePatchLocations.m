function [res, imSize, patchSize] = generatePatchLocations(...
    imSize, patchSize, varargin)
% generatePatchLocations Split image area into patches and return patch
% locations.
%   [res, croppedImSize, patchSize] = generatePatchLocations(imSize, patchSize)
%   generates a matrix of patch locations of size [nPatches, 2], giving the
%   top-left corner of each patch in image coordinates, assuming that the
%   patches are non-overlapping, of size `patchSize` (which can be a scalar
%   or a pair in the order [nRows, nCols]), and that the image
%   size is `imSize` (again, in the order rows, columns). The options below
%   can be used to employ overlapping patches instead. The second return
%   argument indicates the size of the image after its right
%   and bottom edges were potentially trimmed so that an integer number of
%   patches fits within it.
%
%   The function also returns an 'extended' `patchSize`, which is a
%   two-element vector mirroring the input patch size. If the input is a
%   scalar, then the output contains two copies of that scalar. If the
%   input `patchSize` is empty, then the returned value is equal to `imSize`.
%
%   generatePatchLocations(imSize, patchSize, stride) uses the given stride
%   to generate potentially overlapping patches. The stride can be a scalar
%   or a pair, just like `patchSize`.
%
%   Set `patchSize` to an empty matrix to generate a single patch for the
%   entire image.
%
%   Options:
%    'maxPatchesPerImage': integer
%       If given, this ensures that the number of patches per image doesn't
%       exceed the given value. If it would, random sampling without
%       replacement is used to select `maxPatchesPerImage` of them to keep.
%       (default: no maximum)
%
%   The output is a structure with the following fields:
%    'patchLocations': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels.
%    'gridSize': [nRowsInGrid, nColsInGrid]
%       Number of patches along each dimension.
%    'patchSize':
%    'stride':
%    'maxPatchesPerImage':
%       These are just copies of the input arguments and options.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('stride', [], @(v) isempty(v) || (isnumeric(v) && ...
    isvector(v) && length(v) >= 1 && length(v) <= 2));

parser.addParameter('maxPatchesPerImage', inf, @(x) isnumeric(x) && isscalar(x) && x >= 0);

% display defaults
if nargin == 1 && strcmp(imSize, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% preprocess some parameters

% handle square or rectangular patches and strides
if isscalar(patchSize)
    patchSize = [patchSize patchSize];
end
if isscalar(params.stride)
    params.stride = [params.stride params.stride];
elseif isempty(params.stride)
    % empty stride is the same as it being equal to patchSize
    params.stride = patchSize;
end

if ~isempty(patchSize)
    % figure out where the patches are
    
    % crop image so that we get only full patches
    nPatches = 1 + floor((imSize - patchSize) ./ params.stride);
    imSize = params.stride .* (nPatches - 1) + patchSize;
    
    % generate a grid with the given stride
    [locsRowP, locsColP] = ind2sub(nPatches, 1:prod(nPatches));
    locsP = [locsRowP(:) locsColP(:)];
    locs = 1 + bsxfun(@times, locsP - 1, params.stride(:)');
else
    % whole image is a patch
    if ~isempty(params.stride)
        error([mfilename ':badargs'], 'Can''t use a nontrivial stride without a patch size.');
    end
    nPatches = [1 1];
    locs = [1 1];
    patchSize = imSize;
end

% make sure we don't have too many patches
if ~isempty(params.maxPatchesPerImage) && size(locs, 1) > params.maxPatchesPerImage
    subIdxs = randperm(size(locs, 1), params.maxPatchesPerImage);
    locs = locs(subIdxs, :);
end

% fill the output structure
res = struct;
res.patchLocations = locs;
res.patchSize = patchSize;
res.stride = params.stride;
res.maxPatchesPerImage = params.maxPatchesPerImage;
res.gridSize = nPatches;

end