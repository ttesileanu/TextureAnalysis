function [imgPatches, numWide, numHigh, xcoords, ycoords] = patchifyImage(img, patchSize)
% patchifyImage Split image into patches.
%   imgPatches = patchifyImage(img) splits the image 'img' into patches of
%   size 'patchSize', returning a 3d array of size [patchSize, patchSize, nPatches].
%
%   [imgPatches, numWide, numHigh] = patchifyImage(...) also returns the
%   number of patches fitting horizontally and vertically in the image,
%   numWide and numHigh. The number of patches is nPatches = numWide*numHigh.
%
%   [..., xcoords, ycoords] = patchifyImage(...) also returns the location
%   of the patches, in units of patchSize, i.e., [xcoords(i), ycoords(i)]
%   give the x and y positions of the ith patch.

numHigh = floor(size(img, 1) / patchSize);
numWide = floor(size(img, 2) / patchSize);
nPatches = numWide*numHigh;

imgPatches = zeros(patchSize, patchSize, nPatches);

% get rid of incomplete patches
img = img(1:patchSize*numHigh, 1:patchSize*numWide);
imgRes = reshape(img, patchSize, numHigh, patchSize*numWide);
for j = 1:numHigh
    imgPatches(:, :, ((j-1)*numWide+1):(j*numWide)) = ...
        reshape(imgRes(:,j,:), patchSize, patchSize, numWide);
end

xcoords = zeros(nPatches, 1);
ycoords = zeros(nPatches, 1);
for j = 1:numWide
    for k = 1:numHigh
        l = j + (k-1)*numWide; % this is the patch index for postion [j, k]
        xcoords(l) = j;
        ycoords(l) = k;
    end
end

end