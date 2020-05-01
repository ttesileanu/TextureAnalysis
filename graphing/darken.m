function cdark = darken(c, a)
% DARKEN Return a darker color.
%   cdark = DARKEN(c, a) returns a darker version of the input RGB color.
%   This linearly interpolates between `c` and black, with `a = 0` mapping
%   to `c` and `a = 1` mapping to black.

cdark = (1 - a) * c;
cdark = min(max(cdark, 0), 1);

end