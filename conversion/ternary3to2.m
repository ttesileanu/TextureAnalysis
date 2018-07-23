function vecs2 = ternary3to2(vecs3)
% TERNARY3TO2 Convert 3-component (quasi-probability) vectors in ternary
% texture groups to 2-component locations in a space in which the
% probability simplex is equilateral.
%   vecs2 = TERNARY3TO2(vecs3) converts an n x 3 matrix of locations
%   (quasi-probabilities) in ternary texture groups to an n x 2 matrix of
%   locations in the 2d plane in which the simplex is represented by an
%   equilateral triangle.
%
%   This function assumes that each row of `vecs3` sums to 1.

% handle empty inputs
if isempty(vecs3)
    vecs2 = [];
    return;
end

vecs2 = vecs3*[-1/2 -sqrt(3)/2 ; 1 0 ; -1/2 sqrt(3)/2];

end