function img = reconstitute(patches, extents)
% RECONSTITUTE Reconstitute an image from patches.
%   img = RECONSTITUTE(patches, extents) takes a 3d array of `patches` and
%   an nPatches x 4 matrix of patch extents (in the format [row1, col1,
%   row2, col2]), and builds an image from them. Places that are not
%   covered by any patches are filled with zeros, while places where
%   several patches overlap keep only the data from the patch with the
%   highest index.

nRows = max(extents(:, 3));
nCols = max(extents(:, 4));

img = zeros(nRows, nCols);

for i = 1:size(patches, 3)
    crtExt = extents(i, :);    
    img(crtExt(1):crtExt(3), crtExt(2):crtExt(4)) = patches(:, :, i);
end

end