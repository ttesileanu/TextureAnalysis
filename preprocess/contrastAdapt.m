function [imgOut, crop] = contrastAdapt(img, varargin)
% contrastAdapt Do contrast adaptation on an image, globally and/or locally.
%   imgOut = contrastAdapt(img) performs contrast adaptation on the image
%   by z-scoring its intensity values, and then rescaling the result to map
%   [-1, 1] to [0, 1]. A random result (drawn i.i.d. from a normal
%   distribution with mean and standard deviation equal to 0.5) is returned
%   if the standard deviation is exactly zero. If `img` is a 3d array, the
%   contrast adaptation is applied separately to each slice in the 3rd
%   dimension.
%
%   contrastAdapt(..., 'clip', true) clips the final result to [0, 1].
%
%   imgOut = contrastAdapt(img, patchSize) calculates the value for each
%   pixel by z-scoring within the patch of size `patchSize` surrounding it.
%   The patches surrounding pixels close to the edges of the image are
%   cropped to fit within the image.
%
%   NOTE: This latter option is very slow compared to contrast adapting on
%         non-overlapping patches! The slowdown factor is dependent on patch
%         size. You can expect it to be around (pixels per patch)/4, so
%         about 250x slower for a 32x32 patch.
%
%   [imgOut, crop] = contrastAdapt(...) also returns a matrix with four
%   columns, each row of which identifying the range in the original image
%   patches corresponding to the equalized patches. The format of each row
%   is [row1, col1, row2, col2]. A row that is equal to all zeros
%   corresponds to patches that are skipped in `imgOut` (see 'minLevels'
%   below).
%
%   Options:
%    'minStd'
%       If the standard deviation within the patch is smaller than this
%       given number, it  will be ignored (i.e., skipped in `imgOut`, and
%       the corresponding row in `crop` will be all zeros).
%    'meanFct'
%       Function to use to calculate the center value for contrast
%       adaptation (z-scoring corresponds to mean).

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('patchSize', [], @(p) isempty(p) || (isnumeric(p) && isscalar(p) && p >= 1) || ...
    (isnumeric(p) && isvector(p) && numel(p) == 2 && all(p >= 1)));
% parser.addParameter('perpixel', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('clip', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('minStd', 1e-6, @(n) isscalar(n) && isnumeric(n));
parser.addParameter('meanFct', @mean, @(f) isa(f, 'function_handle'));

% defaults
if nargin == 1 && strcmp(img, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

nPatches = size(img, 3);
crop = repmat([1 1 size(img, 1) size(img, 2)], nPatches, 1);
if isempty(params.patchSize)
    % already patchified image
    imgOut = zeros(size(img));
    patchMask = true(nPatches, 1);
    for i = 1:nPatches
        patch = img(:, :, i);
        patchStd = std(patch(:));
        if patchStd < params.minStd
            patchMask(i) = false;
            continue;
        else
            patch = (patch - params.meanFct(patch(:))) / patchStd;
        end
        
        imgOut(:, :, i) = 0.5*(1 + patch);
    end
    imgOut = imgOut(:, :, patchMask);
    crop(~patchMask, :) = 0;
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
                patchStd = std(patch(:));
                if patchStd == 0
                    imgOut(i, j, k) = 0.5*(1 + randn);
                else
                    imgOut(i, j) = 0.5*(1 + (img(i, j) - params.meanFct(patch(:))) / patchStd);
                end
            end
        end
    end
end

if params.clip
    imgOut = min(max(imgOut, 0), 1);
end

end