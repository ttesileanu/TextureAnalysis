function M = fitEllipse(data)
% fitEllipse Fit ellipse or ellipsoid to data.
%   M = fitEllipse(data) returns an n x n matrix `M` identifying the
%   ellipsoid that most closely approximates the `N` `n`-dimensional
%   datapoints in the N x n `data` matrix. The ellipse is given by the
%   equation
%       x'*M*x = 1.
%
%   Note that the objective function is based on how far `x'*M*x` is from
%   1, which is not the same as the Euclidean distance from the points to
%   the ellipsoid.

% figure out in how many dimensions we're working today
[~, n] = size(data);

% find covariance matrix, focus on independent components
rhs = unrollSymmetric(data'*data);

% XXX should document what I'm doing here, but I already forgot...
A = zeros(n*(n+1)/2);
k = 1;
for q = 1:n
    for p = q:n
        if p ~= q
            factor = 2;
        else
            factor = 1;
        end
        A(:, k) = factor*unrollSymmetric(data'*diag(data(:, p) .* data(:, q))*data);
        k = k + 1;
    end 
end

M = rollSymmetric(A\rhs);

end

function v = unrollSymmetric(m)
% Unroll a symmetric matrix into a vector with independent components.

% check the inputs
n = size(m, 1);
if ~ismatrix(m) || size(m, 2) ~= n
    error('This shouldn''t happen: non-symmetric matrix input to unrollSymmetric.');
end

% XXX this is slow and dumb...
v = zeros(n*(n+1)/2, 1);
k = 1;
for j = 1:n
    for i = j:n
        v(k) = m(i, j);
        k = k + 1;
    end
end

end

function m = rollSymmetric(v)
% Roll back a symmetric matrix from its vector representation.

p = length(v);
sqrtDisc = sqrt(1 + 8*p);
n = (-1 + sqrtDisc)/2;
if 1 + 8*p < 0 || ~isvector(v) || abs(floor(n) - n) > eps
    error('This shouldn''t happen: vector doesn''t have the right number of elements to be rolled to a symmetric matrix.');
end

m = zeros(n);
k = 1;
for j = 1:n
    for i = j:n
        m(i, j) = v(k);
        m(j, i) = m(i, j);
        k = k + 1;
    end
end

end