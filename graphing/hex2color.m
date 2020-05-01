function rgb = hex2color(html)
% HEX2COLOR Convert hex color to Matlab-format RGB.
%   rgb = HEX2COLOR(html) converts a color given in web (hex) format into a
%   3-element RGB vector with entries from 0 to 1.

rgb = hex2dec(reshape(html, [], 3)')'/255;

end
