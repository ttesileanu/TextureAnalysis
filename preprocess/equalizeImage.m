function [imgOut, crop] = equalizeImage(img, varargin)
% equalizeImage Perform histogram equalization in an image, potentially
% patch-by-patch.
%   imgOut = equalizeImage(img) performs histogram equalization on the image
%   by mapping all intensity values to the interval [0, 1], and ensuring that
%   they are uniformly spaced. Identical values are separated by adding a
%   small amount of random noise to all pixels (see below to change or
%   deactivate this).
%
%   imgOut = equalizeImage(img, patchSize) splits the image into
%   non-overlapping patches of size `patchSize` and then equalizes within
%   each patch. `patchSize` can also be a pair of integers to use non-square
%   patches. In this case, the order is row (y) size first, and then column
%   (x) size, `[patchSizeY, patchSizeX]`! Anything left over on the right
%   and bottom edges of the image if `size(img)` is not an integer multiple
%   of `patchSize` is trimmed.
%
%   imgOut = equalizeImage(img, patchSize, 'perpixel', true) calculates the
%   value for each pixel by equalizing within the patch of size `patchSize`
%   surrounding it. The patches surrounding pixels close to the edges of
%   the image are cropped to fit within the image.
%
%   NOTE: The 'perpixel' option is slow by a patch-size-dependent factor.
%         You can expect the factor to be of the order
%         (pixels per patch)/25, so about 40x slower for a 32x32 patch.
%
%   [imgOut, crop] = equalizeImage(...) also returns a vector of four
%   elements [row1, col1, row2, col2] identifying the range in the original
%   image corresponding to the equalized image.
%
%   equalizeImage(..., 'jitter', amount) sets the amount of jitter added to
%   the image to distinguish pixels with exactly equal values. The jitter is
%   added as uniform noise in the interval [-amount, amount]. Set this to 0
%   to have no jitter. The default value is equal to `eps(max(abs(img(:))))`.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('patchSize', [], @(p) isempty(p) || (isnumeric(p) && isscalar(p) && p >= 1) || ...
    (isnumeric(p) && isvector(p) && numel(p) == 2 && all(p >= 1)));
parser.addParameter('perpixel', false, @(b) islogical(b) && isscalar(b));
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

% handle defaults
if isempty(params.patchSize)
    if params.perpixel
        error([mfilename ':needpatch'], 'The perpixel option requires a patch size.');
    end
    params.patchSize = size(img);
elseif numel(params.patchSize) == 1
    params.patchSize = [params.patchSize params.patchSize];
end

if isempty(params.jitter)
    params.jitter = eps(max(abs(img(:))));
end

% jitter, if necessary
if params.jitter > 0
    img = img + (2*params.jitter)*(rand(size(img)) - 0.5);
end

if ~params.perpixel
    % non-overlapping patches
    nPatches = floor(size(img) ./ params.patchSize);
    imgOut = zeros(params.patchSize .* nPatches);
    
    for i = 1:nPatches(1)
        ys = 1 + (i-1)*params.patchSize(1);    % starting row for patch
        rows = ys:ys+params.patchSize(1)-1;    % all rows within patch
        for j = 1:nPatches(2)
            xs = 1 + (j-1)*params.patchSize(2); % starting column for patch
            cols = xs:xs+params.patchSize(2)-1; % all columns within patch
            
            patch = img(rows, cols);
            [~, idxs] = sort(patch(:));
            patch(idxs) = linspace(0, 1, length(idxs));
            
            imgOut(rows, cols) = patch;
        end
    end
else
    % a patch around each pixel
    imgOut = zeros(size(img));
    
    dr = floor(params.patchSize/2);    
    [ymax, xmax] = size(img);

    for i = 1:ymax
        % rows in surrounding patch; clip at edges
        rows = max(1, i-dr(1)):min(ymax, i+dr(1));
        for j = 1:xmax
            % columns in surrounding patch; clip at edges
            cols = max(1, j-dr(2)):min(xmax, j+dr(2));
            patch = img(rows, cols);
            
            imgOut(i, j) = sum(patch(:) < img(i, j)) / (numel(patch) - 1);
        end
    end
end

crop = [1 1 size(imgOut)];

end