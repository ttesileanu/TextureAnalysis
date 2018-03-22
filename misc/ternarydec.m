function [mags, uvecs] = ternarydec(vecs)
% TERNARYDEC Decompose a set of vectors in ternary texture planes into
% magnitudes times vectors located on the unit circle.
%   [mags, uvecs] = TERNARYDEC(vecs) decomposes the n x 3 matrix of
%   locations in some ternary texture plane into an n-component vector of
%   magnitudes, `mags`, and an n x 3 matrix of 'normalized' vectors,
%   `uvecs`. The vectors are normalized such that the direction from the
%   [1/3, 1/3, 1/3] origin of the coordinate system to the vector is
%   preserved, and the L2 norm is unity.
%
%   This function assumes that each row of the input vector sums to 1.

% w = v - 1/3
% w0 = a*w
% v0 = w0 + 1/3
%   s.t. v0'*v0 = 1
% (a*(v'-1/3) + 1/3)*(a*(v-1/3) + 1/3) = (a*v' + 1/3*(1-a))*(a*v + 1/3*(1-a))
%       = a^2*(v'*v) + 2/3*(1-a)*a + 1/3*(1-a)^2
% (assuming sum(v) = 1)
%       = a^2*(v'*v) + 2/3*a - 2/3*a^2 + 1/3 - 2/3*a + 1/3*a^2
%       = a^2*((v'*v) - 1/3) + 1/3
%      != 1
%
% So, a^2*((v'*v) - 1/3) = 2/3
%     a^2 = 2/(3*(v'*v) - 1)

vecs0 = vecs - 1/3;
vmags2 = sum(vecs.^2, 2);
mags = sqrt((3*vmags2 - 1) / 2);
uvecs0 = bsxfun(@rdivide, vecs0, mags);
uvecs = 1/3 + uvecs0;

end