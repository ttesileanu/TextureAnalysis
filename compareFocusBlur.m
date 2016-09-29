function [patches, patch_idxs] = compareFocusBlur(m, n, res, varargin)
% compareFocusBlur Make a figure comparing in-focus to blurred patches.
%   compareFocusBlur(m, n, res) generates a figure comparing m*n in-focus
%   image patches to m*n blurred ones. The patches can be chosen at random
%   or deterministically (see the options) from the contents of the `res`
%   structure.
%
%   More specifically, `res.focus.clusterIds` maps each patch to one of
%   several 'components'. The in-focus component is identified by
%   `res.focus.focusCluster`. `res.imgIds` identifies the image index
%   for the patch, which then maps to `res.imageNames` for the file name.
%   The path is given by `res.path` (or overriden by the 'imgPath' option
%   -- see below). The patch is found using the patch positions from
%   `res.patchLocationsOrig`.
%
%   Options:
%    'choice' <s/v>
%       This can be either the string 'random' or a 2 x m*n matrix giving
%       the indices of the patches to consider for the in-focus and blurred
%       images, respectively (in the latter case, `res.focus.clusterIds`
%       and `res.focus.focusCluster` are ignored).
%       (default: 'random')
%    'imgPath' <s>
%       Override the path from which to load the images.
%       (default: `res.path`)
%    'sharpnessMeasure' <c/v>
%       If provided, use a continuous sharpness measure instead of the
%       cluster identifications to find in-focus and blurred patches. This
%       can be a string -- the name of a subfield of `res.focus` to use as
%       a sharpness measure -- or directly a vector of sharpness values for
%       each patch.
%    'clim' <[cmin, cmax]>
%       If provided, color range used for the patches.
%       (default: different scaling for each patch)

% parse the optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('imgPath', '', @(s) ischar(s) && isvector(s));
parser.addParameter('choice', 'random', @(c) (ischar(c) && isvector(c)) || ...
    (isnumeric(c) && isequal(size(c), [2 m*n])));
parser.addParameter('sharpnessMeasure', '', @(c) (ischar(c) && isvector(c)) || ...
    (isnumeric(c) && isvector(c) && isreal(c)));
parser.addParameter('clim', [], @(clim) isvector(clim) && numel(clim) == 2);

parser.parse(varargin{:});
params = parser.Results;

% default image path
if isempty(params.imgPath)
    params.imgPath = res.path;
end

% handle the sharpness measure, if provided
if ~isempty(params.sharpnessMeasure)
    if ischar(params.sharpnessMeasure)
        params.sharpnessMeasure = res.focus.(params.sharpnessMeasure);
    end
end

% decide on a size for the figure
fig_size = [633*n/(2*m) 600];

figure;
set(gcf, 'position', [50 50 fig_size], 'paperposition', [0.25 0.25 fig_size/100]);

% decide which patches to draw
if strcmp(params.choice, 'random')
    if isempty(params.sharpnessMeasure)
        infocus_idxs = flatten(find(res.focus.clusterIds == res.focus.focusCluster))';
        blurred_idxs = flatten(find(res.focus.clusterIds ~= res.focus.focusCluster))';
    else
        % count as "in-focus" patches that are above the 75th percentile in
        % the sharpness measure
        % and "blurry" patches that are below the 25th percentile
        thresh = quantile(params.sharpnessMeasure, [1/4, 3/4]);
        infocus_idxs = flatten(find(params.sharpnessMeasure >= thresh(2)))';
        blurred_idxs = flatten(find(params.sharpnessMeasure <= thresh(1)))';
    end
    
    % XXX there will be an error if there aren't at least m*n entries in
    % either infocus_idxs or blurred_idxs
    
    infocus_sel = randperm(length(infocus_idxs), m*n);
    blurred_sel = randperm(length(blurred_idxs), m*n);
    
    params.choice = [infocus_idxs(infocus_sel) ; blurred_idxs(blurred_sel)];
end

patch_idxs = reshape(params.choice, [2 m n]);
image_names = res.imageNames(res.imgIds(patch_idxs));

patch_locs = res.patchLocationsOrig(patch_idxs, :);

patches = cell(2, m, n);
for k = 1:numel(image_names)
    image = loadLUMImage(fullfile(params.imgPath, image_names{k}));
    crt_patch = image(patch_locs(k, 1):patch_locs(k, 3), patch_locs(k, 2):patch_locs(k, 4));
    patches{k} = crt_patch;
end

if isempty(params.clim)
    opts = {};
else
    opts = {params.clim};
end
for i = 1:m
    for j = 1:n
        crt = (i-1)*n + j;
        
        subplot(2*m, n, crt);
        imagesc(patches{1, i, j}, opts{:});
        colormap('gray');
        hold on;
        if ~isempty(params.sharpnessMeasure)
            text(size(patches{1, i, j}, 2)/2, size(patches{1, i, j}, 1)/2, ...
                num2str(params.sharpnessMeasure(patch_idxs(1, i, j)), '%.2f'), ...
                'color', 'r');
        end
        
        if i == 1
            title('In-focus');
        end
        
        subplot(2*m, n, crt + m*n);
        imagesc(patches{2, i, j}, opts{:});
        colormap('gray');
        if ~isempty(params.sharpnessMeasure)
            text(size(patches{2, i, j}, 2)/2, size(patches{2, i, j}, 1)/2, ...
                num2str(params.sharpnessMeasure(patch_idxs(2, i, j)), '%.2f'), ...
                'color', 'r');
        end
        
        if i == 1
            title('Blurred');
        end
    end
end

end