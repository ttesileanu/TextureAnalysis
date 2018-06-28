function [imgOut, crop] = equalizeImage(img, varargin)
% equalizeImage Perform histogram equalization in an image, potentially
% patch-by-patch.
%   imgOut = equalizeImage(img) performs histogram equalization on the image
%   by mapping all intensity values to the interval [0, 1], and ensuring that
%   they are uniformly spaced. Identical values are separated by adding a
%   small amount of random noise to all pixels (see below to change or
%   deactivate this). If `img` is a 3d array, the equalization is applied
%   separately to each slice in the 3rd dimension.
%
%   imgOut = equalizeImage(img, patchSize) calculates the value for each
%   pixel by equalizing within the patch of size `patchSize` surrounding it.
%   The patches surrounding pixels close to the edges of the image are
%   cropped to fit within the image.
%
%   NOTE: The 'perpixel' option is very slow compared to equalizing on
%         non-overlapping patches! The slowdown factor is dependent on patch
%         size. You can expect it to be around (pixels per patch)/25, so
%         about 40x slower for a 32x32 patch.
%
%   [imgOut, crop] = equalizeImage(...) also returns a vector of four
%   elements [row1, col1, row2, col2] identifying the range in the original
%   image corresponding to the equalized image. If the input image is 2d,
%   this is always [1 1 size(img)]. If it is 3d, the `crop` is set to an
%   empty matrix.
%
%   Options:
%    'jitter'
%       Set the amount of jitter added to the image to distinguish pixels
%       with exactly equal values. The jitter is added as uniform noise in
%       the interval [-amount, amount]. Set this to 0 to have no jitter.
%       The default value is equal to `eps(max(abs(img(:))))`.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('patchSize', [], @(p) isempty(p) || (isnumeric(p) && isscalar(p) && p >= 1) || ...
    (isnumeric(p) && isvector(p) && numel(p) == 2 && all(p >= 1)));
% parser.addParameter('perpixel', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('jitter', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));

% defaults
if nargin == 1 && strcmp(img, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle jitter
if isempty(params.jitter)
    params.jitter = eps(max(abs(img(:))));
end

% jitter image, if necessary
if params.jitter > 0
    img = img + (2*params.jitter)*(rand(size(img)) - 0.5);
end

nPatches = size(img, 3);
if isempty(params.patchSize)
    % non-overlapping patches
    imgOut = zeros(size(img));
    patch = zeros(size(img, 1), size(img, 2));
    sortedValues = linspace(0, 1, numel(patch));
    % equalizing the histogram is the same as labeling each color value by
    % its rank (and we normalize to keep everything in [0, 1])
    for i = 1:nPatches
        [~, idxs] = sort(flatten(img(:, :, i)));
        patch(idxs) = sortedValues;
        imgOut(:, :, i) = patch;
    end
else
    % handle scalar or vector patch sizes
    if numel(params.patchSize) == 1
        params.patchSize = [params.patchSize params.patchSize];
    end
    
    imgOut = zeros(size(img));
    
    % a patch around each pixel, going dr(1) in row direction, dr(2) in
    % column direction
    dr = floor(params.patchSize/2);
    ymax = size(img, 1);
    xmax = size(img, 2);
    
    for k = 1:nPatches
        for i = 1:ymax
            % rows in surrounding patch; clip at edges
            rows = max(1, i-dr(1)):min(ymax, i+dr(1));
            for j = 1:xmax
                % columns in surrounding patch; clip at edges
                cols = max(1, j-dr(2)):min(xmax, j+dr(2));
                patch = img(rows, cols, k);
                
                imgOut(i, j, k) = sum(patch(:) < img(i, j, k)) / (numel(patch) - 1);
            end
        end
    end
end

if nPatches == 1
    crop = [1 1 size(imgOut)];
else
    crop = [];
end

end