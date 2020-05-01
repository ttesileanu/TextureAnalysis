function clight = lighten(c, a)
% LIGHTEN Return a lighter color.
%   clight = lighten(c, a) returns a lighter version of the input RGB
%   color. This linearly interpolates between `c` and white, with `a = 0`
%   mapping to `c` and `a = 1` mapping to white.

clight = (1 - a) * c + a;
clight = min(max(clight, 0), 1);

end