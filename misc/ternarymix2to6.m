function vecs6 = ternarymix2to6(vecs2, dirs)
% TERNARYMIX2TO6 Convert 2-component coordinates in mixed ternary
% statistics groups to the 6-component (quasi-)probability representation.
%   vecs6 = TERNARYMIX2TO6(vecs2, dirs) converts an n x 2 matrix `vecs2` of
%   locations in a mixed ternary group to the 6-component
%   (quasi-)probability representation, the n x 6 matrix `vecs6`. The
%   directions in the two planes are indicated by the 2-component vector
%   `dirs`. In this vector, 0 represents the probability vector [1, 0, 0],
%   1 represents [0, 1, 0], and 2 represents [0, 0, 1].

% first find the quasi-probabilities in the canonical directions (indicated
% by dirs)
qp2 = (2*vecs2 + 1)/3;

% now find the probabilities in the other directions
other_qp2 = (1 - qp2)/2;

% convert to 6-dimensional vectors
vecs6 = [repmat(other_qp2(:, 1), 1, 3) repmat(other_qp2(:, 2), 1, 3)];
vecs6(:, [1 4] + dirs) = qp2;

end