function res = iif(varargin)
% IIF Inline conditional.
%   res = IIF(cond1, fct1, cond2, fct, ..., true, fctN) runs the first
%   function `fctI` for which the corresponding condition value `condI` is
%   `true`, and returns its output.
%
%   Adapted from
%   https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/.

res = varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

end