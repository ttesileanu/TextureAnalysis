function res = analyzePatchGradients(image, patchSize, varargin)
% analyzePatchGradients Calculate statistics derived from various gradients
% for patches of an image.
%   res = analyzePatchGradients(image, patchSize) splits the image into
%   non-overlapping patches of size `patchSize` and calculates the
%   statistics of images obtained by taking various kinds of gradients (see
%   below) and averaging them or calculating their variances.
%
%   The options below can be used to use overlapping patches. `patchSize`
%   can be a single number or a pair of numbers in order to use rectangular
%   patches. In the latter case, the order is row (y) size first, and then
%   column (x) size, `[patchSizeY, patchSizeX]`! When non-overlapping
%   patches are used, the right and bottom edges of the image are ignored
%   if the image size is not an integer multiple of the patch size.
%
%   The function calculates patch averages and standard deviations for the
%   gradient magnitude (as estimated using a Sobel filter) and the
%   Laplacian.
%
%   analyzePatchGradients(image, patchSize, mask) focuses only on the
%   parts of the image contained within the given `mask`. The options below
%   can be used to choose whether only patches fully-contained within the
%   mask are to be used, or whether partial patches are fine.
%
%   Set `patchSize` to an empty matrix to analyze the entire image, or the
%   part of the image contained within the `mask`.
%
%   Options:
%    'maskCrop':
%    'maxPatchesPerImage':
%    'minPatchUsed':
%    'overlapping':
%       See `generatePatchLocations`.
%
%   The output is a structure with the following fields:
%    'patchLocations':
%    'patchLocationsOrig':
%    'pxPerPatch':
%    'patchSize':
%    'overlapping':
%    'minPatchUsed':
%       Inherited from `generatePatchLocations` output.
%
%    'ev': [nPatches, nStats] matrix
%       Matrix containing the calculated statistics for each of the
%       patches.
%    'statsDesc': [nStats, 1] cell array
%       Cell array of strings with human-readable descriptions of the
%       statistics that are returned.

% find patch locations
[res, ~, patchSize, matchedMask] = generatePatchLocations(size(image), ...
    patchSize, varargin{:});

% set up the gradient calculations
gradientFcts = {@getSobelGradients, @getLaplacian};
summaryFcts = {@getSummaryMean, @getSummaryStd};

% calculate gradients on the entire image
% XXX should probably remove edges
gradDesc = {};
imGradients = {};
for k = 1:length(gradientFcts)
    [crtGrad, crtGradDesc] = gradientFcts{k}(image);
    imGradients = [imGradients crtGrad]; %#ok<AGROW>
    gradDesc = [gradDesc crtGradDesc]; %#ok<AGROW>
end

% calculate stats for the gradients
ev = [];
nPatches = size(res.patchLocations, 1);
statsDesc = {};
for i = 1:nPatches
    crtLoc = res.patchLocations(i, :);
    rows = (crtLoc(1):crtLoc(1) + patchSize(1) - 1);
    cols = (crtLoc(2):crtLoc(2) + patchSize(2) - 1);
    maskPatch = matchedMask(rows, cols);
    crtEv = [];
    for j = 1:length(imGradients)
        patch = imGradients{j}(rows, cols);
        % replacing pixels that are not to be used by NaNs, so they are
        % ignored in the analysis
        patch(~maskPatch) = nan;
        for k = 1:length(summaryFcts)
            [crtSummary, crtSummaryDesc] = summaryFcts{k}(patch);
            crtEv = [crtEv crtSummary]; %#ok<AGROW>
            if i == 1
                crtFullDesc = [gradDesc{j} '_' crtSummaryDesc];
                statsDesc = [statsDesc crtFullDesc]; %#ok<AGROW>
            end
        end
    end
    if i == 1
        % allocate space only once, for speed
        ev = zeros(nPatches, numel(crtEv));
    end
    ev(i, :) = crtEv; %#ok<AGROW>
end

% fill the output structure
res.ev = ev;
res.statsDesc = statsDesc;

end

function [grads, names] = getSobelGradients(image)
% Apply Sobel filter, return x, y components and magnitude.

[gx, gy] = imgradientxy(image, 'sobel');
gmag = imgradient(gx, gy);

grads = {gx gy gmag};
names = {'sobel_x', 'sobel_y', 'sobel_mag'};

end

function [grads, names] = getLaplacian(image)
% Apply laplacian filter.

lapFlt = [0 1 0; 1 -4 1; 0 1 0];
filtered = conv2(image, lapFlt, 'same');

grads = {filtered};
names = {'laplace'};

end

function [summary, name] = getSummaryMean(patch)
% Calculate mean ignoring NaNs.

summary = nanmean(patch(:));
name = 'avg';

end

function [summary, name] = getSummaryStd(patch)
% Calculate standard deviation ignoring NaNs.

summary = nanstd(patch(:));
name = 'std';

end