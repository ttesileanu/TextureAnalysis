function [res, imSize, patchSize, matchedMask] = generatePatchLocations(imSize, patchSize, varargin)
% generatePatchLocations Split image area into patches constrained by a
% mask and return patch locations.
%   [res, croppedImSize, patchSize] = generatePatchLocations(imSize, patchSize)
%   generates a matrix of patch locations of size [nPatches, 2], giving the
%   top-left corner of each patch in image coordinates, assuming that the
%   patches are non-overlapping, of size `patchSize` (which can be a scalar
%   or a pair in the order *[patchSizeY, patchSizeX]*), and that the image
%   size is `imSize` (again, in the order *[imSizeY, imSizeX]*). The
%   options below can be used to employ overlapping patches instead. The
%   second return argument indicates the size of the image after its right
%   and bottom edges were potentially trimmed so that an integer number of
%   patches fits within it.
%
%   The function also returns the `patchSize`, which is a two-element
%   vector mirroring the input patch size. If the input is a scalar, then
%   the output contains two copies of that scalar. If the input `patchSize`
%   is empty, then the returned value is equal to `imSize`.
%
%   [..., matchedMask] = generatePatchLocations(imSize, patchSize, mask)
%   focuses only on the parts of the image contained within the given
%   `mask`. The options below can be used to choose whether only patches
%   fully-contained within the mask are to be used, or whether partial
%   patches are fine. The mask can be of a different size than the image
%   (see 'maskCrop' option below). The second return argument is set to a
%   mask that is cropped and resized to match the image size.
%
%   Set `patchSize` to an empty matrix to generate a single patch for the
%   entire image.
%
%   Options:
%    'maskCrop': [row1, col1, row2, col2]
%       Crop region within the mask that corresponds to the image. If
%       [row2 - row1 + 1, col2 - col1 + 1] is not equal to size(image),
%       this can also represent a scaling.
%    'maxPatchesPerImage': integer
%       If given, this ensures that the number of patches per image doesn't
%       exceed the given value. If it would, random sampling without
%       replacement is used to select `maxPatchesPerImage` of them to keep.
%       (default: no maximum)
%    'minPatchUsed': double
%       Minimum fraction of patch that should be contained in the mask.
%       Patches that overlap with the mask less than this fraction are
%       excluded from the analysis. Set this to 1 to only include patches
%       that are fully-contained within the mask..
%       (default: 0)
%    'overlapping': logical, numeric, or pair of numbers
%       If set to true, 1, or (1, 1), the function uses overlapping patches,
%       one for each pixel of the input image (provided it's contained in
%       the `mask`, if one is used). In this case, patches are centered on
%       the pixels. Patches that would go outside the area of the image are
%       not considered (so the pixels close to the edges do not generate
%       patches). This can also be set to an integer larger than 1, or a
%       pair of integers [stepY, stepX], which enables a different spacing
%       to be used between the pixels that generate patches. In all cases
%       the patches are centered on the chosen pixels. A step of 0, or a
%       value of false of this option leads to the use of non-overlapping
%       patches.
%       (default: false)
%
%   The output is a structure with the following fields:
%    'patchLocations': [nPatches, 2] matrix
%       Locations of the top-left corner of the patches, in pixels.
%    'patchLocationsOrig': [nPatches, 4] matrix
%       Locations of the patches ([row1, col1, row2, col2]) in unscaled
%       mask coordinates, if the `maskCrop` option is used (see above).
%       Otherwise this is simply implied by `patchLocations` and
%       `patchSize`.
%    'pxPerPatch': vector
%       The number of pixels that are contained in the `mask` in each patch.
%    'patchSize':
%    'overlapping':
%    'minPatchUsed':
%       These are just copies of the input arguments and options.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('mask', [], @(m) isempty(m) || (ismatrix(m) && isreal(m) && ~ischar(m)));

parser.addParameter('maxPatchesPerImage', inf, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('minPatchUsed', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
parser.addParameter('overlapping', false, ...
    @(b) (isscalar(b) && (islogical(b) || (isnumeric(b) && b >= 0))) || ...
         (isnumeric(b) && isvector(b) && all(b) >= 0 && length(b) == 2));
parser.addParameter('maskCrop', [], @(c) isempty(c) || (isvector(c) && isreal(c) && numel(c) == 4 && all(c >= 1)));

% parse
parser.parse(varargin{:});
params = parser.Results;

% preprocess some parameters
if any(params.overlapping == 0)
    params.overlapping = [];
elseif islogical(params.overlapping) && params.overlapping
    params.overlapping = [1 1];
elseif isnumeric(params.overlapping) && isscalar(params.overlapping)
    params.overlapping = repmat(params.overlapping, 1, 2);
end
params.overlapping = round(params.overlapping);

% default mask is to use all image
if isempty(params.mask)
    params.mask = true(imSize);
    % mask crop is meaningless if there's no mask
    params.maskCrop = [1 1 size(params.mask)];
    maskStep = 1;
else
    % default crop: entire mask
    if isempty(params.maskCrop)
        params.maskCrop = [1 1 size(params.mask)];
        maskStep = 1;
    else
        % figure out scaling from image to mask coordinates
        rowStep = (diff(params.maskCrop([1 3])) + 1)/imSize(1);
        colStep = (diff(params.maskCrop([2 4])) + 1)/imSize(2);
        if abs(rowStep - colStep) >= 1e-6
            error([mfilename ':badmaskcrop'], 'Scaling from mask to image should be aspect-ratio preserving.');
        else
            maskStep = max(1, floor(rowStep));
        end
        
        % scale mask to image coordinates, making sure that only pixels that
        % are fully contained in the mask are counted
        params.mask = params.mask(params.maskCrop(1):params.maskCrop(3), ...
            params.maskCrop(2):params.maskCrop(4));
        params.mask = (blockAverage(double(params.mask), maskStep, 'avg') >= 0.999999);
    end
end
if ~all(imSize == size(params.mask))
    error([mfilename ':badmasksize'], 'Incompatible mask and image sizes.');
end

% handle square or rectangular patches
if numel(patchSize) == 1
    patchSize = [patchSize patchSize];
end

if ~isempty(patchSize)
    % figure out where the patches might be
    if isempty(params.overlapping)
        % ensure image size is integer multiple of patch size
        nPatches = floor(imSize ./ patchSize);
        imSize = [nPatches(1)*patchSize(1), nPatches(2)*patchSize(2)];
        
        [locsRowP, locsColP] = ind2sub(nPatches, 1:prod(nPatches));
        locs = [(locsRowP(:)-1)*patchSize(1) (locsColP(:)-1)*patchSize(2)] + 1;
    else
        % find patch centers -- one for each pixel within the mask,
        % with spacing giving by the steps in params.overlapping
        [locsRow0, locsCol0] = ind2sub(...
            floor(imSize./params.overlapping), ...
            1:prod(imSize));
        % convert from centers to corners
        locsRow = 1 + (locsRow0-1)*params.overlapping(1) - floor(patchSize(1)/2);
        locsCol = 1 + (locsCol0-1)*params.overlapping(2) - floor(patchSize(2)/2);
        % keep only patches that don't go outside the image bounds
        locsMask = ((locsRow >= 1) & (locsCol >= 1) & ...
            (locsRow + patchSize(1) - 1 <= imSize(1)) & ...
            (locsCol + patchSize(2) - 1 <= imSize(2)));
        locsRow = locsRow(locsMask);
        locsCol = locsCol(locsMask);
        locs = [locsRow(:) locsCol(:)];
    end
    wholeImagePatch = false;
else
    % whole image is a patch
    if ~isempty(params.overlapping)
        error([mfilename ':badargs'], 'Can''t use overlapping patches without a patch size.');
    end
    locs = [1 1];
    patchSize = imSize;
    wholeImagePatch = true;
end

% find the patches that have enough overlap with the mask
locsMask = true(size(locs, 1), 1);
pxPerPatch = zeros(size(locs, 1), 1);

locsOrig0 = bsxfun(@plus, params.maskCrop(1:2), (locs-1)*maskStep);
locsOrig = [locsOrig0 bsxfun(@plus, locsOrig0, maskStep*patchSize(:)' - 1)];
for i = 1:size(locs, 1)
    crtLoc = locs(i, :);
    rows = (crtLoc(1):crtLoc(1) + patchSize(1) - 1);
    cols = (crtLoc(2):crtLoc(2) + patchSize(2) - 1);
    maskPatch = params.mask(rows, cols);
    
    % skip patches that don't have enough useable pixels, unless the whole
    % image is a patch
    nPix = sum(maskPatch(:));
    pxPerPatch(i) = nPix;
    if ~wholeImagePatch && nPix/numel(maskPatch) < params.minPatchUsed
        locsMask(i) = false;
    end
end

% restrict to only the patches that were valid
locs = locs(locsMask, :);
locsOrig = locsOrig(locsMask, :);
pxPerPatch = pxPerPatch(locsMask, :);

% make sure we don't have too many patches
if ~isempty(params.maxPatchesPerImage) && imSize(1) > params.maxPatchesPerImage
    subIdxs = randperm(imSize(1), params.maxPatchesPerImage);
    locs = locs(subIdxs, :);
    locsOrig = locsOrig(subIdxs, :);
    pxPerPatch = pxPerPatch(subIdxs, :);
end

% fill the output structure
res = struct;
res.patchLocations = locs;
res.patchLocationsOrig = locsOrig;
res.pxPerPatch = pxPerPatch;
res.patchSize = patchSize;
res.overlapping = params.overlapping;
res.maxPatchesPerImage = params.maxPatchesPerImage;
res.minPatchUsed = params.minPatchUsed;

matchedMask = params.mask;

end