% generate (parts of the) plot describing method

%% Choose an image and a region to focus on

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');
img0 = loadLUMImage(fullfile('NaturalImages', images{1}));
img0log = convertToLog(img0);

img = contrastAdapt(img0log(3*128:5*128-1, 4*128:6*128-1));

figure;
imagesc(img); axis equal; colormap('gray');

%% Load filter

filter0 = open(fullfile('filters', 'filter4x32.mat'));
filter = filter0.filter;

%% Draw

% make a figure of the right size
fig = figure;
fig.Units = 'pixels';
fig.Position = [30 20 900 154];

% make axes filling the entire figure
ax = axes;
ax.Position = [0 0 1 1];
ax.YDir = 'reverse';

hold on;

% draw only the images, we'll add labels and arrows in post-processing
img_brighten = @(m) 0.5*m + 0.4;

% select whether to draw patches with border or without
dimage = @image;
% dimage = @bimage;

% original image
top_edge = 24;
left_edge = 24;
patch_dist = 24;
panel_dist = (fig.Position(3)*2 - 256*5 - 2*left_edge - patch_dist/2) / 4;
% h1 = image(left_edge, top_edge, img/max(img(:)), 'CDataMapping', 'scaled');
h1 = dimage(left_edge, top_edge, img_brighten(img), 'CDataMapping', 'scaled');
colormap('gray');

% down-sampled (block-averaged) image
img_small = blockAverage(img, 4);
left_small = left_edge + size(img, 2) + panel_dist;
h2 = dimage([left_small left_small + size(img, 2)-1], ...
    [top_edge top_edge + size(img, 1)-1], img_brighten(img_small), 'CDataMapping', 'scaled');
% colormap('gray');

% draw patches
left_patches = left_small + size(img, 2) + panel_dist;
for i = 0:1
    if i == 0
        top_patch = top_edge - patch_dist/2;
    else
        top_patch = top_edge + size(img, 1)/2 + patch_dist/2;
    end
    for j = 0:1
        if j == 0
            left_patch = left_patches - patch_dist / 2;
        else
            left_patch = left_patches + size(img, 2)/2 + patch_dist/2;
        end
        crt_patch = img_small(1+i*32:(i+1)*32, 1+j*32:(j+1)*32);
        dimage([left_patch left_patch + size(img, 2)/2-1], ...
            [top_patch top_patch + size(img, 1)/2-1], img_brighten(crt_patch), ...
            'CDataMapping', 'scaled');
    end
end

% draw whitened patches
img_filtered0 = contrastAdapt(filterImage(img_small, filter), size(filter));
% make sure the range for img_filtered is similar to img
img_filtered = img_filtered0;
left_whitened = left_patches + size(img, 2) + panel_dist;
for i = 0:1
    if i == 0
        top_patch = top_edge - patch_dist/2;
    else
        top_patch = top_edge + size(img, 1)/2 + patch_dist/2;
    end
    for j = 0:1
        if j == 0
            left_patch = left_whitened - patch_dist / 2;
        else
            left_patch = left_whitened + size(img, 2)/2 + patch_dist/2;
        end
        crt_patch = img_filtered(1+i*32:(i+1)*32, 1+j*32:(j+1)*32);
        dimage([left_patch left_patch + size(img, 2)/2-1], ...
            [top_patch top_patch + size(img, 1)/2-1], img_brighten(crt_patch), ...
            'CDataMapping', 'scaled');
    end
end

% draw ternarized patches
img_ternarized = quantize(img_filtered0, 3);
left_ternarized = left_whitened + size(img, 2) + panel_dist;
for i = 0:1
    if i == 0
        top_patch = top_edge - patch_dist/2;
    else
        top_patch = top_edge + size(img, 1)/2 + patch_dist/2;
    end
    for j = 0:1
        if j == 0
            left_patch = left_ternarized - patch_dist / 2;
        else
            left_patch = left_ternarized + size(img, 2)/2 + patch_dist/2;
        end
        crt_patch = img_ternarized(1+i*32:(i+1)*32, 1+j*32:(j+1)*32);
        dimage([left_patch left_patch + size(img, 2)/2-1], ...
            [top_patch top_patch + size(img, 1)/2-1], crt_patch, ...
            'CDataMapping', 'scaled');
    end
end

% no need for axis labeling
axis off;

% set axis units to pixels
xlim(ax, [0 2*fig.Position(3)]);
ylim(ax, [0 2*fig.Position(4)]);

% set the clim
ax.CLim = [0 1];

% set the background to white
fig.Color = [1 1 1];

% save
fname = fullfile('figs', 'draft', 'description_patches.png');
fig_frame = getframe(fig);
pixels = fig_frame.cdata;
imwrite(pixels, fname, 'png');

%% Store a single patch

% NOTE: need to run previous cell before this!

crt_patch = img_ternarized(1:32, 1:32);
% imagesc(crt_patch, [0 1]);
bimage(0, 0, crt_patch, 'CDataMapping', 'scaled');
set(gca, 'clim', [0 1]);
colormap('gray');
axis equal;
axis off;

preparegraph;

safe_print(fullfile('figs', 'draft', 'example_patch'));

%% Make and store a histogram and glider ticks

% NOTE: need to run previous cells before this!

% generate the gliders
gliders = cell(1, 81);
crt = 1;
for D = 0:2
    for C = 0:2
        for B = 0:2
            for A = 0:2
                gliders{crt} = [A B ; C D];
                
                crt = crt + 1;
            end
        end
    end
end

% get the stats
crt_histo = processBlock(crt_patch, 3);

% draw the histogram
fig = figure;
fig.Units = 'inches';
% fig.Position = [2 2 2 1];
fig.Position = [2 2 8 4];

ax = axes;
ax.OuterPosition = [0 0.2 1 0.8];

% crt_color = parula(2); crt_color = crt_color(1, :);
crt_color = [0.5 0.7 1]/2;
bar(crt_histo, 0.6, 'edgecolor', 'none', 'facecolor', crt_color);
xlim([0 82]);

box off;
ax.XTick = [];
ax.YTick = [];

ax_pos = ax.Position;

% draw the gliders along the x axis
ax2 = axes;
ax2.Position = [ax.Position(1) 0 ax.Position(3) ax.Position(2)];
ax2.YDir = 'reverse';
hold on;

glider_size = 2;
top_edge = 10;
for i = 1:81
    row = mod(i-1, 3);
    bimage(i - glider_size/2, top_edge + row*glider_size*1.5, gliders{i}, ...
        'CDataMapping', 'scaled', 'gridwidth', 1);
end
colormap('gray');

xlim([0 82]);
% ylim([0 ax2.Position(4)/ax2.Position(3)*82]);

ax2.CLim = [0 2];

axis equal;
axis off;

% for some reason Matlab defaults to 'opengl' for this, which has
% resolution problems
% fig.Renderer = 'painters';

preparegraph;

safe_print(fullfile('figs', 'draft', 'example_histo'));

%% Show ternarization

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.1 0.7];

crt_patch_filtered = img_filtered(1:32, 1:32);
[patch_hist, patch_bins] = hist(crt_patch_filtered(:), 32);
patch_cdf = cumsum(patch_hist) / sum(patch_hist);

hold on;

cutoffs_y = [1/3, 2/3];
cutoffs_x = zeros(1, 2);
cutoffs_idx = zeros(1, 2);
for i = 1:length(cutoffs_y)
    crt_idx = find(patch_cdf > cutoffs_y(i), 1);
    crt_frac = (cutoffs_y(i) - patch_cdf(crt_idx-1)) / (patch_cdf(crt_idx) - patch_cdf(crt_idx-1));
    cutoffs_x(i) = (1-crt_frac)*patch_bins(crt_idx-1) + crt_frac*patch_bins(crt_idx);
    cutoffs_idx(i) = crt_idx - 1;
end

% fill([patch_bins(1) cutoffs_x(1) cutoffs_x(1) patch_bins(cutoffs_idx(1):-1:1)], ...
%      [0 0 cutoffs_y(1) patch_cdf(cutoffs_idx(1):-1:1)], 'k');
% fill([cutoffs_x(1) cutoffs_x(2) cutoffs_x(2) patch_bins(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_x(1)], ...
%      [0 0 cutoffs_y(2) patch_cdf(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_y(1)], ...
%      [0.5 0.5 0.5]);

fill([patch_bins(1) patch_bins(1) cutoffs_x(1) patch_bins(cutoffs_idx(1):-1:1)], ...
     [0 cutoffs_y(1) cutoffs_y(1) patch_cdf(cutoffs_idx(1):-1:1)], 'k', 'edgecolor', 'none');
fill([patch_bins(1) patch_bins(1) cutoffs_x(2) patch_bins(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_x(1)], ...
     [cutoffs_y(1) cutoffs_y(2) cutoffs_y(2) patch_cdf(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_y(1)], [0.5 0.5 0.5], ...
     'edgecolor', 'none');

for i = 1:length(cutoffs_x)
    plot([cutoffs_x(i) patch_bins(end)], repmat(cutoffs_y(i), 1, 2), 'k:');
end

plot(patch_bins, patch_cdf, 'color', [0    0.4470    0.7410]);

xlim(patch_bins([1 end]));
tl = get(gca, 'ticklength');
set(gca, 'xtick', [], 'ytick', [0 1/3 2/3 1], 'yticklabel', {'0', '1/3', '2/3', '1'}, ...
    'yminortick', 'off', 'ticklength', [0 tl(2)]);

xlabel('luminance');
ylabel('cdf');

beautifygraph('fontscale', 0.5, 'ticksize', 12);
preparegraph;

safe_print(fullfile('figs', 'draft', 'ternarization_cdf'));

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.1 0.7];

hold on;
bar(patch_bins(1:cutoffs_idx(1)), patch_hist(1:cutoffs_idx(1)), 'k', 'linewidth', 0.2);
bar(patch_bins(cutoffs_idx(1)+1:cutoffs_idx(2)), patch_hist(cutoffs_idx(1)+1:cutoffs_idx(2)), ...
    'facecolor', [0.5 0.5 0.5], 'linewidth', 0.2);
bar(patch_bins(cutoffs_idx(2)+1:end), patch_hist(cutoffs_idx(2)+1:end), 'w', 'linewidth', 0.2);

xlabel('luminance');
ylabel('count');

tl = get(gca, 'ticklength');
set(gca, 'xtick', [], 'ytick', max(patch_hist)*[0 1/3 2/3 1], 'yticklabel', {'0', '1/3', '2/3', '1'}, ...
    'yminortick', 'off', 'ticklength', [0 tl(2)]);
xlim(patch_bins([1 end]));
ylim([0 max(patch_hist)]);

beautifygraph('fontscale', 0.5, 'ticksize', 12);
preparegraph;

safe_print(fullfile('figs', 'draft', 'ternarization_pdf'));