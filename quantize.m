function [imgOut, crop] = quantize(img, n, patchSize, type)
% quantize Quantize images to discrete color levels.
%   imgOut = quantize(img, n) discretizes the matrix `img` to `n` levels,
%   such that each occurs approximately `numel(img)/n` times.
%
%   imgOut = quantize(img, n, patchSize) splits the image into
%   non-overlapping patches of size `patchSize` and then discretizes within
%   each patch. `patchSize` can also be a pair of integers to use
%   non-square patches. In this case, the order is row (y) size first,
%   and then column (x) size, `[patchSizeY, patchSizeX]`! Anything left over
%   on the right and bottom edges of the image if `size(img)` is not an
%   integer multiple of `patchSize` is trimmed.
%
%   imgOut = quantize(img, n, patchSize, 'perpixel') calculates the value for
%   each pixel by quantizing within the patch of size `patchSize`
%   surrounding it. The patches surrounding pixels close to the edges of
%   the image are cropped to fit within the image.
%
%   NOTE: the 'perpixel' option is very slow!
%
%   [imgOut, crop] = quantize(...) also returns a vector of four elements
%   [row1, col1, row2, col2] identifying the range in the original image
%   corresponding to the filtered image.

if nargin < 3
    type = 'full';
else
    % handle square or rectangular patches
    if numel(patchSize) == 1
        patchSize = [patchSize patchSize];
    end
    
    if nargin < 4
        type = 'patch';
    end
end

crop = [1 1 size(img)];

switch type
    case 'full'
        if n > 2
            thresholds = quantile(img(:), n);
        elseif n == 2
            thresholds = quantile(img(:), 1/2);
        else
            error([mfilename ':badn'], 'Number of levels should be at least 2.');
        end
        
        imgOut = zeros(size(img));
        for level = 1:n-2
            imgOut(img >= thresholds(level) & img < thresholds(level+1)) = level;
        end
        imgOut(img >= thresholds(end)) = n-1;
    case 'patch'
        nPatches = floor(size(img) ./ patchSize);
        imgOut = zeros(patchSize .* nPatches);
        for i = 1:nPatches(1)
            ys = 1 + (i-1)*patchSize(1);    % starting row for patch
            rows = ys:ys+patchSize(1)-1;    % all rows within patch
            for j = 1:nPatches(2)
                xs = 1 + (j-1)*patchSize(2); % starting column for patch
                cols = xs:xs+patchSize(2)-1; % all columns within patch
                imgOut(rows, cols) = quantize(img(rows, cols), n);
            end
        end
        crop(3:4) = size(imgOut);
    case 'perpixel'
        imgOut = zeros(size(img));
        dr = floor(patchSize/2);
        
        [ymax, xmax] = size(img);
        
        % precalculate the quantiles we will use
        p = (1/n:1/n:1);
        for i = 1:ymax
            yr = max(1, i-dr(1)):min(ymax, i+dr(1));
            for j = 1:xmax
                xr = max(1, j-dr(2)):min(xmax, j+dr(2));
                patch = img(yr, xr);
                thresholds = quantile(patch(:), p);
                imgOut(i, j) = find(img(i, j) <= thresholds, 1) - 1;
            end
        end   
    otherwise
        error([mfilename ':badtype'], 'Unrecognized type.');
end

end