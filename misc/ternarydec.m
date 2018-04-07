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
%
%   The function also works in a similar way if the input has 3*nplanes
%   columns: the vectors are normalized around the [1/3, ..., 1/3] origin
%   such that their L2 norm is fixed. In this case, each row must sum to
%   nplanes, and the normalized norm is sqrt((2 + nplanes) / 3). That is,
%   the norm is chosen such that a vector that has norm 1 in one plane and
%   is at the origin in all the others will be kept unchanged.

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

if mod(size(vecs, 2), 3) ~= 0
    error([mfilename ':badsz'], 'The input argument must have a number of columns that is divisible by 3.');
end

nplanes = size(vecs, 2) / 3;

vecs0 = vecs - 1/3;
vmags2 = sum(vecs.^2, 2);
mags = sqrt((3*vmags2 - nplanes) / 2);
uvecs0 = bsxfun(@rdivide, vecs0, mags);
uvecs = 1/3 + uvecs0;

end