function m = chop_ev(m, threshold)
% CHOP_EV Chop matrix eigenvalues to a threshold while presevering
% eigenvectors.
%   M = CHOP_EV(m0) return a matrix with the same eigenvectors as `m0`, but
%   replacing all negative eigenvalues with 0.
%
%   M = CHOP_EV(m0, threshold) clips all eigenvalues smaller than
%   `threshold` to that threshold.

if nargin < 2
    threshold = 0;
end

[V, D] = eig(m);

D(D < threshold) = threshold;

m = V*D/V;

end