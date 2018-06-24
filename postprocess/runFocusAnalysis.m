function resOut = runFocusAnalysis(res, varargin)
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
%    'skipSharpness': boolean
%       Set to true to skip sharpness calculation. If this is true, then a
%       `focusPatch` or `focusImage` must be provided.
%    'quiet': boolean
%       Set to true to inhibit any display from the function.
%
%   See also: getSharpnessMeasure.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('focusImage', [], @(n) isscalar(n) && isnumeric(n) && isreal(n) && n >= 1);
parser.addParameter('focusPatch', [], @(n) isscalar(n) && isnumeric(n) && isreal(n) && n >= 1);
parser.addParameter('skipSharpness', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('quiet', false, @(b) isscalar(b) && islogical(b));

% defaults
if nargin == 1 && strcmp(res, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% get sharpness for each patch -- need to load the images
t0 = tic;
if ~params.skipSharpness
    allImgIds = unique(res.imgIds);
    focus = struct;
    focus.sharpness = zeros(size(res.ev, 1), 1);
    % cumbersome to use walkImages here, so doing it manually...
    if ~params.quiet
        progress = TextProgress('checking patch sharpness');
    end
    for i = 1:length(allImgIds)
        crtImgId = allImgIds(i);
        crtImg = loadLUMImage(res.imageNames{crtImgId});
        crtPatchMask = (res.imgIds == crtImgId);
        crtPatchLocs = res.patchLocationsOrig(crtPatchMask, :);
        crtSharpness = zeros(size(crtPatchLocs, 1), 1);
        for j = 1:size(crtPatchLocs, 1)
            crtPatchLoc = crtPatchLocs(j, :);
            crtPatch = crtImg(crtPatchLoc(1):crtPatchLoc(3), crtPatchLoc(2):crtPatchLoc(4));
            crtSharpness(j) = getSharpnessMeasure(crtPatch);
        end
        focus.sharpness(crtPatchMask) = crtSharpness;
        if ~params.quiet
            progress.update(100*i/length(allImgIds));
        end
    end
    if ~params.quiet
        progress.done('done!');
    end
else
    if isempty(params.focusImage) && isempty(params.focusPatch)
        error([mfilename ':needidx'], 'When skipSharpness is true, either focusImage or focusPatch must be provided.');
    end
end

% run 2 Gaussian decomposition
if ~params.quiet
    disp('Fitting 2-Gaussian mixture, this might take a while...');
end
options = statset('Display', 'off');
focus.gMix = gmdistribution.fit(res.ev, 2, ...
    'options', options, 'SharedCov', false, 'replicates', 10);

% assign every image patch to one of the two gaussians
focus.clusterIds = focus.gMix.cluster(res.ev);
focus.clusterDistances = focus.gMix.mahal(res.ev);

% for every image, compute average of cluster IDs over all patches
focus.imageClusters = zeros(length(res.imageNames), 1);
for i = 1:length(res.imageNames)
    crtMask = (res.imgIds == i);
    focus.imageClusters(i) = mean(focus.clusterIds(crtMask));
end

% identify the in-focus cluster
if ~isempty(params.focusPatch)
    % easy: was given an in-focus patch
    focus.focusCluster = focus.clusterIds(params.focusPatch);
elseif ~isempty(params.focusImage)
    % for in-focus image, take the cluster that most patches in the image
    % identify with
    focus.focusCluster = round(focus.imageClusters(params.focusImage));
else
    % find the cluster that has the highest median sharpness value
    sharpClusters = arrayfun(@(id) median(focus.sharpness(focus.clusterIds == id)), ...
        1:max(focus.clusterIds));
    [~, focus.focusCluster] = max(sharpClusters);
end

if ~params.quiet
    disp(['Focus analysis took ' num2str(toc(t0), '%.1f') ' seconds.']);
end

% store the focus structure in the output
resOut = res;
resOut.focus = focus;

end