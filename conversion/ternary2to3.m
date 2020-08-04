function vecs3 = ternary2to3(vecs2)
% TERNARY2TO3 Convert 2-component locations in a space in which the
% probability simplex is equilateral to 3-component (quasi-probability)
% vectors in ternary texture groups.
%   vecs3 = TERNARY2TO3(vecs2) converts an n x 2 matrix of locations in the
%   2d plane in which the simplex is represented by an equilateral triangle
%   to an n x 3 matrix of quasi-probabilities in ternary texture groups.
%
%   Each row of `vecs2` will sum to 1.

% handle empty inputs
if isempty(vecs2)
    vecs3 = [];
    return;
end

% u = (2 * y - x - z) / 2
% v = (z - x) * sqrt(3) / 2

% x + y + z = 1 --> y = 1 - x - z

% u = (2 - 2 * x - 2 * z - x - z) / 2 = 1 - 3 * (x + z) / 2
% z + x = 2 * (1 - u) / 3
% z - x = 2 * v / sqrt(3)

% z = (1 - u) / 3 + v / sqrt(3) = 1 / 3 - u / 3 + v / sqrt(3)
% x = (1 - u) / 3 - v / sqrt(3) = 1 / 3 - u / 3 - v / sqrt(3)
% y = 1 - (x + z) = 1 - 2 * (1 - u) / 3 = 1 / 3 + 2 * u / 3

vecs3 = 1/3 + vecs2 * [-1/3 2/3 -1/3 ; -1/sqrt(3) 0 1/sqrt(3)];

end