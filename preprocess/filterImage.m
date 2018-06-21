function [image, crop] = filterImage(image, filter, method)
% filterImage Apply a filter to the image, either on a patch-by-patch basis
% or to the full image at once.
%   filtered = filterImage(image, filter) applies the given filter to the
%   image. The filter must be given as a matrix, and will be applied by
%   first splitting the image into non-overlapping patches of size
%   `size(filter)`, and then filtering each patch. The right and bottom
%   edges of the image are clipped if the image size is not an integer
%   multiple of the filter size.
%
%   filterImage(image, filter, 'full') applies the filter to the whole
%   image instead. An edge of size equal to half the size of the filter is
%   clipped on all four sides of the image to remove edge effects.
%
%   [filtered, crop] = filterImage(...) also returns a vector of four
%   elements [row1, col1, row2, col2] identifying the range in the original
%   image corresponding to the filtered image.

if nargin < 3
    method = 'patch';
end

crop = [1 1 size(image)];

switch method
    case 'patch'
        patchifier = ImagePatchifier({size(image)}, size(filter));
        nPatches = patchifier.gridSize;
        while patchifier.next
            crtCoords = patchifier.getPatchCoordinates;
            rows = crtCoords(1):crtCoords(3);
            cols = crtCoords(2):crtCoords(4);
            image(rows, cols) = real(ifft2(fft2(image(rows, cols)) .* filter));
        end
        
        % XXX the manual version is about 10% faster...
%         nPatches = floor(size(image) ./ size(filter));
%         for i = 1:nPatches(1)
%             ys = 1 + (i-1)*size(filter, 1);    % patch starts here
%             rows = (ys:ys+size(filter, 1)-1);  % all patch rows
%             for j = 1:nPatches(2)
%                 xs = 1 + (j-1)*size(filter, 2);   % patch starts here
%                 cols = (xs:xs+size(filter, 2)-1); % all patch columns
%                 image(rows, cols) = real(ifft2(fft2(image(rows, cols)) .* filter));
%             end
%         end
        
        % trim overflow
        crop(3:4) = nPatches .* size(filter);
        image = image(1:crop(3), 1:crop(4));
    case 'full'
        trimEdge = floor(size(filter)/2);
        cShift = trimEdge - 1;
        image = imfilter(image, circshift(ifft2(filter), cShift));
        
        crop(1:2) = trimEdge;
        crop(3:4) = size(image) - trimEdge + 1;
        image = image(crop(1):crop(3), crop(2):crop(4));
    otherwise
        error([mfilename ':baderror'], 'Unrecognized filtering method.');
end

end