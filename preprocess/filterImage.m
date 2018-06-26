function [image, crop] = filterImage(image, filter)
% filterImage Apply a filter to a 2d or 3d (patchified) array.
%   filtered = filterImage(image, filter) applies the given filter to the
%   image. If the `filter` matrix is the same size as the first two
%   dimensions of `image`, it is applied to each one of the patches
%   identified by the third dimension. If `image` is 2d and has a different
%   size than the `filter`, then after the filte is applied, an edge of
%   size equal to half the size of the filter is clipped on all four sides
%   of the image to remove edge effects.
%
%   [filtered, crop] = filterImage(...) also returns a vector of four
%   elements [row1, col1, row2, col2] identifying the range in the original
%   image corresponding to the filtered image, in the case in which the
%   input image was 2d.

imSize = size(image);

if length(imSize) == 2
    crop = [1 1 imSize];
    nPatches = 1;
else
    crop = [];
    nPatches = imSize(3);
end

if isequal(imSize(1:2), size(filter))
    for i = 1:nPatches
        image(:, :, i) = real(ifft2(fft2(image(:, :, i)) .* filter));
    end
else
    if length(imSize) ~= 2
        error([mfilename ':badsz'], 'Patchifed images require the filter size to match the patch size.');
    end
    
    trimEdge = floor(size(filter)/2);
    cShift = trimEdge - 1;
    image = imfilter(image, circshift(ifft2(filter), cShift));
    
    crop(1:2) = trimEdge;
    crop(3:4) = imSize - trimEdge + 1;
    image = image(crop(1):crop(3), crop(2):crop(4));
end

end