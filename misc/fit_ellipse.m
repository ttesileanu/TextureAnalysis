function M = fit_ellipse(data)
% FIT_ELLIPSE Fit ellipse to data.
%   M = FIT_ELLIPSE(data) returns an nxn matrix `M` identifying the
%   ellipsoid that most closely approximates the N n-dimensional datapoints
%   in the N x n `data` matrix. The ellipse is given by the equation
%       x'*M*x = 1.

[~, n] = size(data);

rhs = unroll_symmetric(data'*data);

A = zeros(n*(n+1)/2);
k = 1;
for q = 1:n
    for p = q:n
        if p ~= q
            factor = 2;
        else
            factor = 1;
        end
        A(:, k) = factor*unroll_symmetric(data'*diag(data(:, p) .* data(:, q))*data);
        k = k + 1;
    end 
end

M = roll_symmetric(A\rhs);

end

function v = unroll_symmetric(m)
% Unroll a symmetric matrix into a vector with independent components.

n = size(m, 1);
if ~ismatrix(m) || size(m, 2) ~= n
    error('This shouldn''t happen: non-symmetric matrix input to unroll_symmetric.');
end

v = zeros(n*(n+1)/2, 1);
k = 1;
for j = 1:n
    for i = j:n
        v(k) = m(i, j);
        k = k + 1;
    end
end

end

function m = roll_symmetric(v)
% Rollback a symmetric matrix from its vector representation.

p = length(v);
sqrt_disc = sqrt(1 + 8*p);
n = (-1 + sqrt_disc)/2;
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