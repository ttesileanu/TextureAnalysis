function subRes = selectRes(res, mask)
% selectRes Select a subset of patches in a results structure.
%   subRes = selectRes(res, mask) returns a structure containing the
%   statistics for a subset of patches in `res`, as identified by `mask`.
%   The mask can be logical or a vector of integers.

subRes = cropSelected(res, {'objIds', 'ev', 'pxPerPatch', ...
    'patchLocations', 'patchLocationsOrig', 'imgIds'}, mask);
if isfield(subRes, 'focus')
    subRes.focus = cropSelected(subRes.focus, ...
        {'sharpness', 'clusterIds', 'clusterDistances'}, mask);
end

end

function s = cropSelected(s, fields, mask)
% Crop selected fields according to the mask. Cropping is done along first
% dimension for matrices and higher-dimensional arrays.

for i = 1:length(fields)
    field = fields{i};
    crtValue = s.(field);
    if isvector(s.(field))
        s.(field) = crtValue(mask);
    else
        oldSize = size(s.(field));
        if islogical(mask)
            newHeight = sum(mask);
        else
            newHeight = length(mask);
        end
        s.(field) = reshape(crtValue(mask, :), [newHeight oldSize(2:end)]);
    end
end

end