function image = applyNonlinearity(image, map)
% applyNonlinearity Apply a nonlinear filter to each pixel of an image.
%   result = applyNonlinearity(image, map) applies a nonlinearity to each
%   pixel of an image. The values in the input image are assumed to be in
%   the range [0, 1] (all others will be clipped), and the `map` is a
%   vector of values corresponding to the nonlinearity. If `f(x)` is the
%   function applying the nonlinearity, then
%       map(i) = f( (i - 1) / (length(map) - 1) ) .
%   In other words, `map` gives the values of the nonlinearity at
%   equally-spaced points in the interval [0, 1], with the ends included.
%   Linear interpolation is used for values that are in-between the ones
%   given in the `map`.

% pad input to make the interpolation code simpler
n_levels = length(map);
map = [map(:) ; map(end)];

% clip input to [0, 1]
image = max(min(image, 1), 0);

% scale to positions in map
image = image*(n_levels-1) + 1;

% interpolate
image_int = floor(image);
image_map0 = map(image_int);
image_map1 = map(image_int + 1);

image_frac = image - image_int;
image = (1 - image_frac) .* image_map0 + image_frac .* image_map1;

end