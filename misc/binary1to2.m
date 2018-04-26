function vecs2 = binary1to2(vecs1)
% BINARY1TO2 Convert 1-component locations in binary texture group space
% to 2-component (quasi-probability) vectors.
%   vecs2 = BINARY1TO2(vecs1) converts an n x 1 matrix of binary texture
%   coordinates into an n x 2 matrix of 2d locations (quasi-probabilities).
%   The 1-component coordinates are simply the differences between
%   the components of the 2-component vectors.
%
%   The function also works for multi-column inputs, in which case each
%   column is transformed separately, and the results are concatenated.

% x = p(1) - p(0) = 1 - 2*p(0) --> p(0) = (1 - x)/2, p(1) = (1+x)/2.

vecs2 = zeros(size(vecs1, 1), 2*size(vecs1, 2));
for i = 1:size(vecs1, 2)
    vecs2(:, 2*i-1:2*i) = 1/2 + vecs1(:, i)*[-1/2 1/2];
end

end