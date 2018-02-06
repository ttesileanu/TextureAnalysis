function continuousMaps = mapInterpolate(discreteValues, discreteKeys, order)
% mapInterpolate Use polynomial fits to interpolate and extrapolate the
% mean of a stochastic mapping.
%   continuousMaps = mapInterpolate(discreteValues, discreteKeys, order)
%   works on a cell array of cell arrays of matrices, `discreteValues`, and
%   a cell array of key vectors, `discreteKeys`. Each entry
%   `discreteValues{i}` gives the values of a stochastic mapping at some
%   number `m` of locations given by `discreteKeys{i}`. The function then
%   focuses on only the means along rows `mean(discreteValues{i}{j}, 1)`,
%   and generates polynomial fits of order `order` to approximate the
%   dependence of these values on the key. In other words,
%       mean(discreteValues{i}{j}, 1) == function_i of discreteKeys{j} .
%
%   The function returns a cell array of structures, each of which
%   containing the following fields:
%    'function'
%       A function mapping keys to (mean) values. This will be a polynomial.
%    'coefficients'
%       A table of coefficients for the polynomial function. Each row in
%       the table corresponds to an output dimension, and each column
%       corresponds to the coefficients of the polynomial (in `polyfit`
%       order).
%   
%   See also: polyfit.

if length(discreteValues) ~= length(discreteKeys) || ~isvector(discreteValues) || ...
        ~isvector(discreteKeys) || ~iscell(discreteKeys) || ~iscell(discreteValues)
    error([mfilename ':badinp'], 'discreteValues and discreteKeys should be cell vectors of the same length.');
end

continuousMaps = cell(size(discreteValues));
for i = 1:length(discreteValues)
    crtMeans = cell2mat(...
        cellfun(@(m) mean(m, 1), discreteValues{i}', 'uniform', false));
    crtKeys = discreteKeys{i};
    crtCoeffs = zeros(size(crtMeans, 2), order+1);
    for j = 1:size(crtMeans, 2)
        crtCoeffs(j, :) = polyfit(crtKeys(:), crtMeans(:, j), order);
    end
    
    crtMap = struct;
    crtMap.coefficients = crtCoeffs;
    crtMap.function = @(x) cell2mat(arrayfun(@(p) x(:).^p, order:-1:0, 'uniform', false)) * ...
        crtCoeffs';
    
    continuousMaps{i} = crtMap;
end

end