function [vecs2, dirs] = ternary6tomix2(vecs6)
% TERNARY6TOMIX2 Convert 6-component (quasi-)probability representation of
% mixed ternary statistics groups to the 2-dimensional representation.
%   [vecs2, dirs] = TERNARY6TOMIX2(vecs6) converts an n x 6 matrix `vecs6`
%   of (quasi-)probability locations in a mixed ternary group to the
%   2-component representation, the n x 2 matrix `vecs2`. This assumes that
%   the departure from iid texture is only along one of the 3 probability
%   axes in each plane. These axes are inferred by the function and
%   returned in the 2-component vector `dirs`. In this vector, 0 represents
%   the probability vector [1, 0, 0], 1 represents [0, 1, 0], and 2
%   represents [0, 0, 1].
%
%   This is the inverse operation to ternarymix2to6.
%
%   The function also generalizes to a group of `k` planes, where `vecs6`
%   has `3*k` columns and vecs2 has `k` columns. The functioning is similar
%   to the above, and is the inverse of ternarymix2to6 in this case also.
%
%   See also: ternarymix2to6.

% handle empty inputs
if isempty(vecs6)
    vecs2 = [];
    dirs = [];
    return;
end

% we can do this more generally, with the argument being 3k dimensional
if mod(size(vecs6, 2), 3) ~= 0
    error([mfilename ':badsz'], 'Input argument should have a number of columns that is a multiple of 3.');
end
k = size(vecs6, 2)/3;

% start building the output
dirs = zeros(1, k);
vecs2 = zeros(size(vecs6, 1), k);
tol = 1e-8;
for i = 1:k
    crt_vecs = vecs6(:, 3*i-2:3*i);
    
    % first to find dirs: we need to know which component is the odd one out
    crt_diffs = abs(diff([crt_vecs crt_vecs(:, 1)], [], 2)) > tol;
    
    % some directions might be neutral in some coordinate planes, so all
    % the differences for those directions will be 0
    crt_non0 = any(crt_diffs, 2);
    
    % hack: replace all zeros by the first non-zero row, if it exists;
    % otherwise set dirs to nan
    non0_row = find(crt_non0, 1);
    if isempty(non0_row)
        dirs(i) = nan;
    else
        for j = 1:3
            crt_diffs(~crt_non0, j) = crt_diffs(non0_row, j);
        end
        
        if max(abs(std(crt_diffs, [], 1))) > tol
            error([mfilename ':nonunif'], 'All input rows should have the same structure.');
        end
        
        % if crt_diffs(:, 1) == false, then crt_vecs(:, 1) == crt_vecs(:, 2),
        % so special direction is the 3rd component, i.e., dir = 2;
        % crt_diffs(:, 1) == false --> 2
        % crt_diffs(:, 2) == false --> 0
        % crt_diffs(:, 3) == false --> 1
        dirs(i) = mod(find(~crt_diffs(1, :), 1) + 1, 3);
    end
    
    vecs2(:, i) = (3/2)*crt_vecs(:, dirs(i) + 1) - 1/2;
end

end