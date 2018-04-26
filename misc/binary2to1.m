function vecs1 = binary2to1(vecs2)
% BINARY2TO1 Convert 2-component (quasi-probability) vectors in binary
% texture group space to 1-component locations.
%   vecs1 = BINARY2TO1(vecs2) converts an n x 2 matrix of 2d locations
%   (quasi-probabilities) into an n x 1 matrix of binary texture
%   coordinates. The 1-component coordinates are simply the differences
%   between the components of the 2-component vectors.
%
%   The function also works for n x 2k input matrices, in which case each
%   pair of columns is transformed separately, and the results are
%   concatenated.

vecs1 = zeros(size(vecs2, 1), size(vecs2, 2)/2);
for i = 1:size(vecs2, 2)/2
    vecs1(:, i) = vecs2(:, 2*i) - vecs2(:, 2*i-1);
end

end