function img = applyImageColorTransformation(img, G, x, y)
% applyImageColorTransformation Apply an affine color transformation to a
% grayscale image.
%   imgOut = applyImageColorTransformation(img, G, x, y) applies an affine
%   transformation to the image `img` with `G` grayscale levels. The
%   grayscale levels are assumed to be uniformly distributed in the
%   interval [0, 1], including the endpoints. These values are rescaled to
%   the range 0, ..., G-1, after which they are transformed according to
%       v -> x*v + y (mod G)
%   and scaled back to the interval [0, 1].
%
%   The function works whether `img` is a matrix or a 3d array of patches.

img = mod(y + x*round((G-1)*img), G)/(G-1);

end