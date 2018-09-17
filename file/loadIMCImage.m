function m = loadIMCImage(name)
% loadIMCImage Load image in van Hateren IML or IMC format.
%   m = loadIMCImage(name) loads an image from van Hateren's database. This
%   can be either IML or IMC format.
%
%   The resolution is fixed at 1536x1024. The image is automatically
%   transposed so that the height goes along the matrix srows.

f = fopen(name, 'rb', 'ieee-be');
w = 1536;
h = 1024;

m = fread(f, [w h], 'uint16');
m = m';

end