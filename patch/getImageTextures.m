function res = getImageTextures(image, nLevels)
% getImageTextures Calculate texture statistics for an image or a set of 
% image patches.
%   res = getImageTextures(image, nLevels) takes in a 3d array of patches,
%   `image`, and calculates the texture statistics with `nLevels` levels
%   for each patch. The image can also be a 2d array, in which case the
%   whole image is treated as a single patch.
%
%   The output is a structure with the following fields:
%    'ev': [nPatches, nStats] matrix
%       Matrix containing the statistics for each of the patches.
%    'nLevels':
%       Copy of input argument.

% process the patches
ev = [];
nPatches = size(image, 3);
for i = 1:nPatches
    patch = image(:, :, i);
    crtEv = analyzeTexture(patch, nLevels);

    if isempty(ev)
        % allocate space only once, for speed
        ev = zeros(nPatches, numel(crtEv));
    end
    ev(i, :) = crtEv; %#ok<AGROW>
end

% fill the output structure
res.ev = ev;
res.nLevels = nLevels;

end