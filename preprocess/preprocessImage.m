function [image, varargout] = preprocessImage(image0, varargin)
% preprocessImage Preprocess the image by going to log space, block
%   averaging, filtering, and quantizing colors.
%   image = preprocessImage(image0, blockAF) preprocesses the image by
%   taking the logarithm of the entries and block averaging using the given
%   factor. `image0` can either be the filename for an image, or a matrix
%   containing the image data itself. If it is a file name, it is read
%   using `loadLUMImage`.
%
%   preprocessImage(image0, blockAF, 'filter', filter) also performs a
%   filtering with the given `filter` after the logarithm and the block
%   averaging. See the options below for more control over the filtering.
%
%   The function can also perform histogram equalization (see 'equalize'
%   option below), quantization (see 'quantize' option), apply a per-pixel
%   non-linear filter (see 'nonlinearity' option), and can perform the
%   first two of these functions in a per-patch fashion (see 'patchSize'
%   option).
%
%   [image, mask1, mask2, ..., maskN, origImage] = ...
%       preprocessImage(image0, mask1, mask2, ..., maskN, ...) preprocesses
%   the image together with a number of 'masks' -- images that will be
%   downsampled without averaging and cropped if necessary, but not
%   filtered, equalized, or quantized. This is useful to ensure that the
%   masks are aligned with the preprocessed image. Note that since the
%   downsampling is done without any averaging, masks with high-frequency
%   components can lead to confusing results. In that case, it might make
%   more sense to not downsample the masks, but instead use the `crops`
%   output argument (see below) to convert between coordinates in the masks
%   and coordinates in the preprocessed image.
%
%   [..., images] = preprocessImage(...) also returns copies of the image
%   along the preprocessing pipeline in a structure with fields 'original',
%   'log', 'averaged', 'filtered', 'equalized'. If a procedure was not
%   performed, the image is identical to that from the previous step in the
%   pipeline (e.g., if equalization is disabled, then `images.equalized` is
%   just a lazy copy of `images.filtered`).
%
%   [..., crops] = preprocessImage(...) also returns a structure indicating
%   where the images returned along the pipeline fit on the original image
%   coordinates. So, crops.averaged, crops.filtered, and crops.final are
%   vectors of four elements [row1, col1, row2, col2] indicating the ranges
%   of rows and columns in the original image that correspond to each of
%   the steps along the preprocessing pipeline. These can be used to convert
%   between coordinates in these images, and coordinates in the original
%   image. (Note that equalization and quantization use the same crop, so
%   that crops.final refers to either one of them or both).
%
%   Options:
%    'averageType': char
%       This is passed to `blockAverage` to set the type of averaging that
%       is performed.
%       (default: use `blockAverage`s default)
%    'doLog': logical
%       Set to `false` to skip the conversion of the image to logarithmic
%       space.
%       (default: true)
%    'equalize': logical
%       if `true`, the image is histogram-equalized after the filtering
%       (and before a potential quantization).
%       (default: true)
%    'equalizeType': char
%       This is passed to `equalize` to set the type of histogram
%       equalization that is used.
%       (default: use the `equalize` default)
%    'filter': [], or matrix
%       Whitening filter to use after log and block-averaging, but before
%       equalization and/or quantization. If empty, no whitening is
%       performed.
%    'filterType': char
%       This is passed to `filterImage` to set the type of filtering that
%       is performed.
%       (default: use `filterImage`s default)
%    'nonlinearity': [], or vector
%       If non-empty, this instructs preprocessImage to perform a per-pixel
%       nonlinear filtering of the image. This happens after equalization
%       but before quantization. See `applyNonlinearity`.
%       (default: [])
%    'patchSize': [], int, or [int, int]
%       Patch size to use for `quantize` and/or `equalize`. See the
%       documentation for those two functions for details. If empty, the
%       whole image is processed at once (but see 'quantType' below).
%       (default: size(filter) if provided, else [])
%    'quantize': int
%       If non-empty, perform color quantization to the given number of
%       levels after equalization (if any). The quantization algorithm
%       assumes values in the range [0, 1]; everything outside that range
%       is clipped.
%    'quantizeType': char
%       This can be 'deterministic' (in which case `quantize` is used), or
%       'stochastic' (in which case `stochasticBinarize` is used). Note
%       that 'stochastic' can only be used with binary images, so
%       'quantize' must be 2 in this case.
%    'threshold': double
%       Threshold used by the conversion to logarithmic space to avoid
%       taking the logarithm of zero or negative numbers.
%       (default: use `convertToLog`s default)
%
%   See also: loadLUMImage, convertToLog, blockAverage, filterImage, equalize, quantize, stochasticBinarize.

% check for masks
firstNonMaskIdx = find(~cellfun(@(m) isreal(m) && ismatrix(m) && ~isscalar(m), ...
    varargin), 1);
if isempty(firstNonMaskIdx) || firstNonMaskIdx == 1
    masks = {};
else
    masks = varargin(1:firstNonMaskIdx-1);
    varargin = varargin(firstNonMaskIdx:end);
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addRequired('blockAF', @(n) isscalar(n) && isreal(n) && n >= 1);

parser.addParameter('filter', [], @(m) isempty(m) || (ismatrix(m) && isreal(m)));
parser.addParameter('quantize', [], @(n) isempty(n) || (isnumeric(n) && isscalar(n) && n > 1));
parser.addParameter('patchSize', [], @(p) isempty(p) || (isnumeric(p) && isscalar(p) && p >= 1) || ...
    (isnumeric(p) && isvector(p) && numel(p) == 2 && all(p >= 1)));

parser.addParameter('averageType', [], @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('doLog', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('filterType', [], @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('equalize', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('equalizeType', [], @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParameter('nonlinearity', [], @(v) isempty(v) || (isvector(v) && isnumeric(v) && isreal(v)));
parser.addParameter('quantizeType', [], @(s) isempty(s) || (ischar(s) && isvector(s) && ...
    ismember(s, {'deterministic', 'stochastic'})));
parser.addParameter('threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle defaults
if isempty(params.quantizeType)
    params.quantizeType = 'deterministic';
end

% check options
if ~isempty(params.quantize)
    if strcmp(params.quantizeType, 'stochastic') && params.quantize ~= 2
        error([mfilename ':badstochq'], 'Stochastic quantization only works with 2 quantization levels.');
    end
end

% load the image from file, if needed
if ischar(image0)
    image = loadLUMImage(image0);
else
    image = image0;
end

% keep track of the original image
origImage = image;

% take the log
if params.doLog
    if isempty(params.threshold)
        image = convertToLog(image);
    else
        image = convertToLog(image, params.threshold);
    end
end
logImage = image;

% block average
if params.blockAF > 1
    if isempty(params.averageType)
        avgOpts = {};
    else
        avgOpts = {params.averageType};
    end
    image = blockAverage(image, params.blockAF, avgOpts{:});
    
    % block average the masks, if any
    for i = 1:numel(masks)
        if ~isempty(masks{i})
            masks{i} = blockAverage(masks{i}, params.blockAF, 'sub'); %#ok<AGROW>
        end
    end
end
averagedImage = image;

% whitening filter
if ~isempty(params.filter)
    % filter
    if isempty(params.filterType)
        filterOpts = {};
    else
        filterOpts = {params.filterType};
    end
    [image, filterCrop] = filterImage(image, params.filter, filterOpts{:});
    
    % crop the masks, if necessary
    for i = 1:numel(masks)
        if ~isempty(masks{i})
            masks{i} = masks{i}(filterCrop(1):filterCrop(3), filterCrop(2):filterCrop(4)); %#ok<AGROW>
        end
    end
else
    filterCrop = [1 1 size(image)];
end
filteredImage = image;

% histogram equalization
if params.equalize
    if isempty(params.patchSize)
        if ~isempty(params.equalizeType)
            error([mfilename ':badeqtype'], 'Non-default equalization types require a patch size.');
        end
        [image, eqCrop] = equalize(image);
    else        
        if isempty(params.equalizeType)
            quantOpts = {};
        else
            quantOpts = {params.equalizeType};
        end
        [image, eqCrop] = equalize(image, params.patchSize, quantOpts{:});
    end
    
    % crop the masks, if necessary
    for i = 1:numel(masks)
        if ~isempty(masks{i})
            masks{i} = masks{i}(eqCrop(1):eqCrop(3), eqCrop(2):eqCrop(4)); %#ok<AGROW>
        end
    end
else
    eqCrop = [1 1 size(image)];
end
equalizedImage = image;

% apply per-pixel nonlinearity
if ~isempty(params.nonlinearity)
    image = applyNonlinearity(image, params.nonlinearity);
end
nonlinearImage = image;

% color quantization
if ~isempty(params.quantize) && isfinite(params.quantize)
    switch params.quantizeType
        case 'deterministic'
            image = quantize(image, params.quantize);
        case 'stochastic'
            image = stochasticBinarize(image);
    end
end

images.origImage = origImage;
images.logImage = logImage;
images.averagedImage = averagedImage;
images.filteredImage = filteredImage;
images.equalizedImage = equalizedImage;
images.nonlinearImage = nonlinearImage;

% update the crops
crops.averaged = [1 1 size(averagedImage)*params.blockAF];
% filtered crop is in scaled coordinates
crops.filtered = [(filterCrop(1:2)-1)*params.blockAF+1 filterCrop(3:4)*params.blockAF];
% combine filter with quantization crop to obtain coordinates in scaled but
% non-trimmed image
qfCropOrigin = filterCrop(1:2) + eqCrop(1:2) - 1;
qfCrop = [qfCropOrigin qfCropOrigin+eqCrop(3:4)-1];
% scale back to coordinates in original image
crops.final = [(qfCrop(1:2)-1)*params.blockAF+1 qfCrop(3:4)*params.blockAF];

% return the image with the masks, if there are any masks
varargout = [masks {images crops}];

end