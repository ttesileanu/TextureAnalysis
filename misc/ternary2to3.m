function vecs3 = ternary2to3(vecs2)
% TERNARY2TO3 Convert 2-component locations in ternary texture group space
% where the probability simplex is equilateral to 3-component
% (quasi-probability) vectors.
%   vecs3 = TERNARY2TO3(vecs2) converts an n x 2 matrix of locations in the
%   2d plane in which the simplex is represented by an equilateral triangle
%   to an n x 3 matrix of 3d locations (quasi-probabilities) in ternary
%   texture groups.

vecs3 = 1/3 + vecs2*[-1/3 2/3 -1/3 ; -sqrt(3)/3 0 sqrt(3)/3];

end