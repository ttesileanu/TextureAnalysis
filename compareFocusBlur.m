function [patches, patch_idxs] = compareFocusBlur(m, n, res, varargin)
% compareFocusBlur Make a figure comparing in-focus to blurred patches.
%   compareFocusBlur(m, n, res) generates a figure comparing m*n in-focus
%   image patches to m*n blurred ones. The patches can be chosen at random
%   or deterministically (see the options) from the contents of the 'res'
%   structure.
%
%   More specifically, res.cx maps each patch to one of several 'components'.
%   The in-focus component is identified by res.focus.component.
%   res.ic.image identifies the image index for the patch, which then maps
%   to res.ic.name for the file name. The patch is found using the patch
%   positions res.ic.x and res.ic.y, which are expressed in units of
%   res.R*res.N.
%
%   Options:
%    'choice' <s/v>
%       This can be either the string 'random' or a 2 x mn matrix giving
%       the indices of the patches to consider for the in-focus and blurred
%       images, respectively (in the latter case, res.cx and
%       res.focus.component are ignored).
%       (default: 'random')
%    'imgpath' <s>
%       The path from which to load the images.
%       (default: current directory)
%    'clim' <[cmin, cmax]>
%       If provided, color range used for the patches.
%       (default: different scaling for each patch)

% parse the optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('imgpath', '', @(s) ischar(s) && isvector(s));
parser.addParamValue('choice', 'random', @(c) (ischar(c) && isvector(c)) || ...
    (isnumeric(c) && isequal(size(c), [2 m*n])));
parser.addParamValue('clim', [], @(clim) isvector(clim) && numel(clim) == 2);

parser.parse(varargin{:});
params = parser.Results;

% decide on a size for the figure
figsize = [633*n/(2*m) 600];

%figure;
set(gcf, 'position', [50 50 figsize], 'paperposition', [0.25 0.25 figsize/100]);

% decide which patches to draw
if strcmp(params.choice, 'random')
    flatten = @(v) v(:);
    
    infocus_idxs = flatten(find(res.cx == res.focus.component))'; %#ok<FNDSB>
    blurred_idxs = flatten(find(res.cx ~= res.focus.component))'; %#ok<FNDSB>
    
    infocus_sel = randperm(length(infocus_idxs), m*n);
    blurred_sel = randperm(length(blurred_idxs), m*n);
    
    params.choice = [infocus_idxs(infocus_sel) ; blurred_idxs(blurred_sel)];
end

patch_idxs = reshape(params.choice, [2 m n]);
image_names = res.ic.name(res.ic.image(patch_idxs));

patch_size = res.R*res.N;
xs = (res.ic.x(patch_idxs)-1)*patch_size;
ys = (res.ic.y(patch_idxs)-1)*patch_size;

patches = zeros(patch_size, patch_size, 2, m, n);
for k = 1:numel(image_names)
    image_struc = open(fullfile(params.imgpath, image_names{k}));
    image = image_struc.LUM_Image;
    crt_patch = image((ys(k)+1):(ys(k)+patch_size), (xs(k)+1):(xs(k)+patch_size));
    patches(:, :, k) = crt_patch;
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
        imagesc(patches(:, :, 1, i, j), opts{:});
        colormap('gray');
        
        if i == 1
            title('In-focus');
        end
        
        subplot(2*m, n, crt + m*n);
        imagesc(patches(:, :, 2, i, j), opts{:});
        colormap('gray');
        
        if i == 1
            title('Blurred');
        end
    end
end

end