function vflat = flatten(v)
% FLATTEN Flatten an array.
%   vflat = FLATTEN(v) flattens the given array, having the same effect as
%   v(:). The advantage of this function is that it can be applied after an
%   indexing operation, for example one can say flatten(v(:, 1)), while
%   v(:, 1)(:) is not allowed in Matlab.

vflat = v(:);

end