function elem = selectMod(v, i)
% selectMod Select elements from a list, wrapping around the end if the
% requested index is larger than the length.
%   elem = selectMod(v, i) selects the `i`th element from the array `v`. If
%   `i` is larger than `numel(v)`, the index is wrapped around, so that
%   `v(numel(v) + 1) == v(1)`, etc. This works for both cell array and
%   numeric arrays.

% wrap
n = numel(v);
i = 1 + mod(i-1, n);

% select
if iscell(v)
    elem = v{i};
else
    elem = v(i);
end

end