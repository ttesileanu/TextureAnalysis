function setdefault(varname, default)
% SETDEFAULT Set variable in base workspace to a default value, unless it's
% already set.
%   SETDEFAULT(varname, default) checks whether the variable with the given
%   name exists in the base workspace. If it does, the function does
%   nothing. Otherwise, it sets the variable to the `default` value.

% make sure the input is a valid variable name
if ~isvarname(varname)
    error([mfilename ':badvar'], 'The first argument needs to be a valid variable name.');
end

b = evalin('base', ['exist(''' varname ''', ''var'')']);
if ~b
    assignin('base', varname, default);
end

end