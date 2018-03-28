function v = getnth(n, fct, varargin)
% GETNTH Get the nth output of the function, or the nth element of an array.
%   v = GETNTH(n, fct, arg1, ... argN) gets the nth output of the function
%   call fct(arg1, ..., argN).
%
%   v = GETNTH(n, vec), where vec is not a function, gets the nth value in
%   the array vec. This is useful in cases where vec is obtained from an
%   expression that cannot be followed by indexing.
%
% (adapted from here
% http://stackoverflow.com/questions/3710466/how-do-i-get-the-second-return-value-from-a-function-without-using-temporary-var)


if isa(fct, 'function_handle')
    [values{1:n}] = fct(varargin{:});
    v = values{n};
else
    v = fct(n);
end

end