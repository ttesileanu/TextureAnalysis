function res = runFocusAnalysis(res, varargin)
% runFocusAnalysis Fit mixture of Gaussians to image patch statistics and
% attempt to find in-focus patches.
%   resOut = runFocusAnalysis(res) fits a mixture of Gaussians to the image
%   patches from `res`. A sharpness measure based on the Laplacian is also
%   calculated for each patch. The Gaussian corresponding to in-focus
%   patches is identified based on this sharpness measure, unless the
%   'focusImage' or 'focusPatch' options are used (see below).
%
%   The results from the focus analysis are stored in a field `focus` of
%   the output structure `resOut`, with the following fields:
%    'gMix':
%       Matlab Gaussian mixture object describing the fit.
%    'sharpness':
%       Vector of sharpness values for each patch. See `getSharpnessMeasure`.
%    'clusterIds':
%       Vector of integers identifying the Gaussian cluster to which each
%       patch is most likely to belong.
%    'clusterDistances':
%       Matrix in which each column gives the Mahalanobis distance from
%       the patches to each Gaussian.
%    'imageClusters':
%       Vector of mean clusterIds for all patches within each image.
%    'focusCluster':
%       An integer identifying which of the Gaussian clusters is more
%       likely to correspond to in-focus patches.
%
%   Options:
%    'focusImage': integer
%       The index of an image that has most of its patches in focus. If
%       this is provided, it is used (instead of the sharpness measure) to
%       identify which of the Gaussian clusters is most likely to contain
%       the in-focus patches.
%    'focusPatch': integer
%       The index of a patch that is in focus. If this is provided, it is
%       used (instead of the sharpness measure) to identify which of the
%       Gaussian clusters is most likely to contain the in-focus patches.
%
%   See also: getSharpnessMeasure.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('focusImage', [], @(n) isscalar(n) && isnumeric(n) && isreal(n) && n >= 1);
parser.addParameter('focusPatch', [], @(n) isscalar(n) && isnumeric(n) && isreal(n) && n >= 1);

% parse
parser.parse(varargin{:});
params = parser.Results;

% get sharpness for each patch -- need to load the images
allImgIds = unique(res.imgIds);
focus = struct;
focus.sharpness = zeros(size(res.ev, 1), 1);
for i = 1:length(allImgIds)
    crtImgId = allImgIds(i);
    crtImgName = res.imageNames{crtImgId};
    crtImg = loadLUMImage(fullfile(res.path, crtImgName));
    crtPatchMask = (res.imgIds == crtImgId);
    crtPatchLocs = res.patchLocationsOrig(crtPatchMask, :);
    crtSharpness = zeros(size(crtPatchLocs, 1), 1);
    for j = 1:size(crtPatchLocs, 1)
        crtPatchLoc = crtPatchLocs(j, :);
        crtPatch = crtImg(crtPatchLoc(1):crtPatchLoc(3), crtPatchLoc(2):crtPatchLoc(4));
        crtSharpness(j) = getSharpnessMeasure(crtPatch);
    end
    focus.sharpness(crtPatchMask) = crtSharpness;
end

% run 2 Gaussian decomposition on all analyses

options = statset('Display', 'off');
focus.gMix = gmdistribution.fit(res.ev, 2, ...
    'options', options, 'SharedCov', false, 'replicates', 10);

% assign every image patch to one of the two gaussians
focus.clusterIds = focus.gMix.cluster(res.ev);
focus.clusterDistances = focus.gMix.mahal(res.ev);

% for every image, compute how many patches are in cluster 1 or cluster 2
% on average
focus.imageClusters = zeros(length(res.imageNames), 1);
for i = 1:length(res.imageNames);
    crtMask = (res.imgIds == i);
    focus.imageClusters(j) = mean(focus.clusterIds(crtMask));
end

% identify the in-focus cluster
if ~isempty(params.focusPatch)
    focus.focusCluster = focus.clusterIds(params.focusPatch);
elseif ~isempty(params.focusImage)
    focus.focusCluster = round(focus.imageClusters(params.focusImage));
else
    % find the cluster that has the highest median sharpness value
    sharpClusters = arrayfun(@(id) median(focus.sharpness(focus.clusterIds == id)), ...
        1:max(focus.clusterIds));
    [~, focus.focusCluster] = max(sharpClusters);
end

res.focus = focus;

end