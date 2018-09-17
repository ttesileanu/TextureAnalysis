function colorXyz = rgbToXyz(color)
% rgbToXyz Convert gamma-compressed sRGB colors to XYZ.
%   colorXyz = rgbToXyz(color) converts each row of the gamma-compressed
%   sRGB color matrix `color` to XYZ. The components of `color` should be
%   between 0 and 1.

% convert to linearized sRGB
threshold = 0.04045;
mask = (color <= threshold);

normalizer = 12.92;
color(mask) = color(mask) / normalizer;

a = 0.055;
gamma = 2.4;
color(~mask) = ((color(~mask) + a) / (1 + a)) .^ gamma;

% linear transformation to XYZ
m = [0.4124 0.3576 0.1805 ; 0.2126 0.7152 0.0722 ; 0.0193 0.1192 0.9505];
colorXyz = color*m';

end
