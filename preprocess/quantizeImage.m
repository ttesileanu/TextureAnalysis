function [img, crop] = quantizeImage(img, n)
% quantizeImage Quantize images to discrete color levels.
%   imgOut = quantizeImage(img, n) discretizes the matrix `img` to `n`
%   levels, equally spaced between 0 and 1. Values outside the [0, 1] range
%   are clipped to 0 or 1.
%
%   [imgOut, crop] = quantizeImage(...) returns a cropping area, always
%   equal to [1 1 size(img)] if `img` is a matrix, and [] otherwise.

img = min(max(floor(n*img)/(n-1), 0), 1);
if ismatrix(img)
    crop = [1 1 size(img)];
else
    crop = [];
end

end