function [mags, uvecs] = binarydec(vecs)
% BINARYDEC Decompose a set of vectors in binary texture groups into
% magnitudes times vectors located on the unit circle.
%   [mags, uvecs] = BINARYDEC(vecs) converts an n x 2 matrix of locations
%   in a binary group into an n-component vector of magnitudes, `mags`, and
%   an n x 2 matrix of normalized (quasi-)probability vectors, `uvecs`:
%   each row of `uvecs` sums up to 1. The vectors are normalized such that
%   the direction from the [1/2, 1/2] origin of the coordinate system to
%   the vector is preserved, and the L2 norm is unity.
%
%   The function also works in a similar way if the input has 2*groups
%   columns: the vectors are normalized around the [1/2, ..., 1/2] origin
%   such that their L2 norm is fixed. In this case, each row must sum to
%   ngroups, and the normalized norm is sqrt((1 + ngroups) / 2). That is,
%   the norm is chosen such that a vector that has norm 1 in one group and
%   is at the origin in all the others will be kept unchanged.

% w = v - 1/2
% w0 = a*w
% v0 = w0 + 1/2
%   s.t. v0'*v0 = 1
% (a*(v'-1/2) + 1/2)*(a*(v-1/2) + 1/2) = (a*v' + 1/2*(1-a))*(a*v + 1/2*(1-a))
%       = a^2*(v'*v) + (1-a)*a + 1/2*(1-a)^2
% (assuming sum(v) = 1)
%       = a^2*(v'*v) + a - a^2 + 1/2 - a + 1/2*a^2
%       = a^2*((v'*v) - 1/2) + 1/2
%      != 1
%
% So, a^2*((v'*v) - 1/2) = 1/2
%     a^2 = 1/(2*(v'*v) - 1)

% ngroups > 1
% v0'*v0 = a^2*(v'*v) + (1-a)*a*ngroups + 1/2*(1-a)^2*ngroups
% (assuming sum(v) = ngroups)
%       = a^2*(v'*v) + ngroups*(a - a^2 + 1/2 - a + 1/2*a^2)
%       = a^2*(v'*v) + ngroups*(1 - a^2)/2
%       = a^2*((v'*v) - ngroups/2) + ngroups/2
%      != (1 + ngroups) / 2
%
% So, a^2*((v'*v) - ngroups/2) = 1/2
%     a^2 = 1/(2*(v'*v) - ngroups)

if mod(size(vecs, 2), 2) ~= 0
    error([mfilename ':badsz'], 'The input argument must have an even number of columns.');
end

ngroups = size(vecs, 2) / 2;

vecs0 = vecs - 1/2;
vmags2 = sum(vecs.^2, 2);
mags = sqrt(2*vmags2 - ngroups);
uvecs0 = bsxfun(@rdivide, vecs0, mags);
uvecs = 1/2 + uvecs0;

end