function subResults = sampleResults(results, mask)
% sampleResults Sample a subset of patches in a results structure.
%   subResults = sampleResults(results, mask) returns a structure
%   containing the statistics for a subset of patches in `results`, as
%   identified by `mask`. The mask can be logical or a vector of integers.

subResults = sampleSelected(results, {'ev', 'patchLocations', 'imageIds'}, mask);
if isfield(subResults, 'focus')
    subResults.focus = sampleSelected(subResults.focus, ...
        {'sharpness', 'clusterIds', 'clusterDistances'}, mask);
end

end

function s = sampleSelected(s, fields, mask)
% Sample selected fields according to the mask. Sampling is done along
% the first dimension for matrices and higher-dimensional arrays.

for i = 1:length(fields)
    field = fields{i};
    crtValue = s.(field);
    if isvector(crtValue)
        s.(field) = crtValue(mask);
    else
        s.(field) = crtValue(mask, :);
    end
end

end