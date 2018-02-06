function [imgOut, crop] = contrastAdapt(img, varargin)
% contrastAdapt Do contrast adaptation on an image, globally and/or locally.
%   imgOut = contrastAdapt(img) performs contrast adaptation on the image
%   by z-scoring its intensity values. A random result (drawn i.i.d. from
%   a normal distribution with mean and standard deviation equal to 0.5) is
%   returned if the standard deviation is exactly zero.
%
%   imgOut = contrastAdapt(img, patchSize) splits the image into
%   non-overlapping patches of size `patchSize` and then performs z-scoring
%   within each patch. `patchSize` can also be a pair of integers to use
%   non-square patches. In this case, the order is row (y) size first, and
%   then column (x) size, `[patchSizeY, patchSizeX]`! Anything left over on
%   the right and bottom edges of the image if `size(img)` is not an
%   integer multiple of `patchSize` is trimmed.
%
%   contrastAdapt(..., 'clip', true) clips the final result to [0, 1].
%
%   imgOut = contrastAdapt(img, patchSize, 'perpixel', true) calculates the
%   value for each pixel by z-scoring within the patch of size `patchSize`
%   surrounding it. The patches surrounding pixels close to the edges of
%   the image are cropped to fit within the image.
%
%   NOTE: The 'perpixel' option is slow by a patch-size-dependent factor.
%         You can expect the factor to be of the order
%         (pixels per patch)/4, so about 1000x slower for a 64x64 patch.
%
%   [imgOut, crop] = contrastAdapt(...) also returns a vector of four
%   elements [row1, col1, row2, col2] identifying the range in the original
%   image corresponding to the equalized image.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('patchSize', [], @(p) isempty(p) || (isnumeric(p) && isscalar(p) && p >= 1) || ...
    (isnumeric(p) && isvector(p) && numel(p) == 2 && all(p >= 1)));
parser.addParameter('perpixel', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('clip', false, @(b) islogical(b) && isscalar(b));

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
            patch_std = std(patch(:));
            if patch_std == 0
                patch = randn(size(patch));
            else
                patch = (patch - mean(patch(:))) / patch_std;
            end
                        
            imgOut(rows, cols) = 0.5*(1 + patch);
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
            patch_std = std(patch(:));
            if patch_std == 0
                imgOut(i, j) = 0.5*(1 + randn);
            else
                imgOut(i, j) = 0.5*(1 + (img(i, j) - mean(patch(:))) / patch_std);
            end
        end
    end
end

if params.clip
    imgOut = min(max(imgOut, 0), 1);
end

crop = [1 1 size(imgOut)];

end