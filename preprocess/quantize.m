function img = quantize(img, n)
% QUANTIZE Quantize images to discrete color levels.
%   imgOut = QUANTIZE(img, n) discretizes the matrix `img` to `n` levels,
%   equally spaced between 0 and 1. Values outside the [0, 1] range are
%   clipped to 0 or 1.

img = min(max(floor(n*img)/(n-1), 0), 1);

end