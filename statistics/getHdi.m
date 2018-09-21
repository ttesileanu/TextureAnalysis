function [lo, hi] = getHdi(v, p)
% getHdi Get highest density interval from vector.
%   [lo, hi] = getHdi(v) finds the interval `[lo, hi]` that contains at
%   least a fraction `p` of all the points in `v` and has the smallest
%   extent (`hi - lo` is as low as possible). This automatically implies
%   that `lo` and `hi` will coincide with points in the vector.

% number of elements that we must have in the interval
n = ceil(p*length(v));

% sorting the data allows us to search linearly for the interval
v = sort(v);

% find extents of interval depending on where it starts
extents = arrayfun(@(i) v(i+n-1) - v(i), 1:length(v)-n+1);

% find the location for the minimal extent
[~, i0] = min(extents);

lo = v(i0);
hi = v(i0 + n - 1);

end
