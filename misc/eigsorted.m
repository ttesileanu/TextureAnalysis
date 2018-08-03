function [V, D] = eigsorted(X, k, varargin)
% EIGSORTED Get eigenvalues and eigenvectors, sorted in decreasing order of
% eigenvalues.
%   E = EIGSORTED(X) returns the eigenvalues of the matrix in descending
%   order.
%
%   E = EIGSORTED(X, k) returns the top k eigenvalues of X.
%
%   [V, D] = EIGSORTED(X) returns the eigenvectors (the columns of V) and
%   a diagonal matrix of eigenvalues of the matrix X, in decreasing order
%   of the eigenvalues.
%
%   [V, D] = EIGSORTED(X, k) returns only the top k eigenvalues and
%   eigenvectors.
%
%   See also: EIG.

% if no number is given, output all eigenvalues/eigenvectors
if nargin < 2
    k = size(X, 1);
end

if nargout < 2
    % only need eigenvalues
    V = eig(X);
    V = sort(V, 'descend');
    V = V(1:k);
else
    % get eigenvalues and eigenvectors
    [V, D] = eig(X);
    [~, indices] = sort(diag(D), 'descend');
    subindices = indices(1:k);
    V = V(:, subindices);
    D = D(subindices, subindices);

%     for n = 1:k
%         [~, idx] = max(abs(V(:, n)));
%         % this can only be zero if V(:, n) is 0, which shouldn't be
%         % possible... and if it is, multiplying it by 0 won't
%         % change anything anyway
%         factor = sign(V(idx(1), n));
%         V(:, n) = factor*V(:, n);
%     end
end

end