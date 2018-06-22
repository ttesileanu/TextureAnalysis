function res = getImageTextures(image, nLevels, patchSize, varargin)
% getImageTextures Calculate texture statistics for patches of an image.
%   res = getImageTextures(image, nLevels, patchSize) splits the image into
%   non-overlapping patches of size `patchSize` and calculates the texture
%   statistics with `nLevels` levels for each patch. The options below can
%   be used to choose overlapping patches instead. `patchSize` can be a
%   single number or a pair of numbers in order to use rectangular patches.
%   In the latter case, the order is row (y) size first, and then column
%   (x) size, `[patchSizeY, patchSizeX]`! When non-overlapping patches are
%   used, the right and bottom edges of the image are ignored if the image
%   size is not an integer multiple of the patch size. Set `patchSize` to
%   an empty matrix to treat the whole image as a single patch.
%
%   getImageTextures(image, nLevels, patchSize, 'mask', mask) focuses only
%   on the parts of the image where `mask ~= 0`. The options below can be
%   used to choose whether only patches fully-contained within the mask are
%   to be used, or whether partial patches are fine.
%
%   Options:
%    'maskCrop':
%    'maxPatchesPerImage':
%    'minPatchUsed':
%    'stride':
%       See `generatePatchLocations`.
%
%   The output is a structure with the following fields:
%    'patchLocations':
%    'patchLocationsOrig':
%    'pxPerPatch':
%    'patchSize':
%    'stride':
%    'minPatchUsed':
%       Inherited from `generatePatchLocations` output.
%
%    'ev': [nPatches, nStats] matrix
%       Matrix containing the statistics for each of the patches.
%    'nLevels':
%       Copy of input argument.
%
%   See also: generatePatchLocations.

% defaults
if nargin == 1 && strcmp(image, 'defaults')
    disp('All arguments are forwarded to generatePatchLocations.');
    return;
end

% find patch locations
[res, ~, patchSize, matchedMask] = generatePatchLocations(size(image), ...
    patchSize, varargin{:});

% process the patches that have enough overlap with the mask
ev = [];
for i = 1:size(res.patchLocations, 1)
    crtLoc = res.patchLocations(i, :);
    rows = (crtLoc(1):crtLoc(1) + patchSize(1) - 1);
    cols = (crtLoc(2):crtLoc(2) + patchSize(2) - 1);
    patch = image(rows, cols);
    
    % use the mask if it was given
    if ~isempty(matchedMask)
        maskPatch = matchedMask(rows, cols);
    
        % replacing pixels that are not to be used by NaNs, so they are
        % ignored in the analysis
        % XXX does this work for nLevels > 2?
        patch(~maskPatch) = nan;
    end
    
    crtEv = analyzeTexture(patch, nLevels);

    if isempty(ev)
        % allocate space only once, for speed
        ev = zeros(size(res.patchLocations, 1), numel(crtEv));
    end
    ev(i, :) = crtEv; %#ok<AGROW>
end

% fill the output structure
res.ev = ev;
res.nLevels = nLevels;

end