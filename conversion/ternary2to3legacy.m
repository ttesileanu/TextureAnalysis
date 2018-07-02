function vecs3 = ternary2to3legacy(vecs2)
% TERNARY2TO3LEGACY Convert legacy 2-component locations in ternary
% texture groups to the 3-component (probability) representation.
%   vecs3 = TERNARY2TO3LEGACY(vecs2) converts an n x 2 matrix of locations
%   in ternary texture groups to an n x 3 matrix of quasi-probabilities
%   (i.e., each row summing up to 1, but potentially having negative
%   elements). This starts from the legacy 2-component encoding used in the
%   raw psychophysics data, where [0 0] maps to the center of the
%   probability simplex (i.e., all probabilities equal to 1/3), [1 0] maps
%   to [1 0 0], [0 1] maps to [0 1 0], and [-1 -1] maps to [0 0 1]. The
%   transformation is affine.

vecs3 = 1/3 + vecs2*[[2/3 ; -1/3] [-1/3 ; 2/3] [-1/3 ; -1/3]];

end