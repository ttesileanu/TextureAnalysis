function [patched, locations] = patchify(image, patchSize, varargin)
% PATCHIFY Split image into patches.
%   patched = PATCHIFY(image, patchSize) splits the `image` into
%   non-overlapping patches of size `patchSize`. This can be a scalar for
%   square patches, or a pair of numbers for rectangular patches (in the
%   worder [rowSize, colSize]). The patches are returned in the form of a
%   3d array.
%
%   patched = PATCHIFY(image, patchSize, stride) uses the given stride to
%   generate patches that can be overlapping. The `stride` can again be
%   either a scalar or a pair of numbers.
%
%   [patched, locations] = PATCHIFY(...) returns the extents of the
%   patches as an nPatches x 4 array in which each row has the form
%   [row1, col1, row2, col2].
%
%   Other options are passed directly to `generatePatchLocations`. Note
%   that unlike that function, PATCHIFY always returns the patch locations
%   in `image` coordinates, even if the patches are selected according to a
%   mask.
%
%   See also: generatePatchLocations.

[patchStruct, ~, patchSize] = generatePatchLocations(size(image), patchSize, varargin{:});

nPatches = size(patchStruct.patchLocations, 1);
patched = zeros([patchSize nPatches]);

locations = [patchStruct.patchLocations ...
    bsxfun(@plus, patchStruct.patchLocations, patchSize(:)'-1)];

for i = 1:nPatches
    crtLoc = locations(i, :);
    patched(:, :, i) = image(crtLoc(1):crtLoc(3), crtLoc(2):crtLoc(4))  ;
end

end