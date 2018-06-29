% make plots that will be used to explain the ternary texture analysis

%% Save independent texture patches

% this will be used on the tips and in the center of a ternary alloy plot

% keep things reproducible
rng(3434);

% choose the plane, locations in each plane (with labels for saving), and patch size
planes = {'AB_1_2'};
texAxes = {[0 1 0], [0 0 1], [1 0 0], [1/3, 1/3, 1/3]};
locLabels = {'1', '2', '0', 'c'};
patchSize = 128;

% draw everything
for i = 1:length(planes)
    plane = planes{i};
    
    for k = 1:length(texAxes)
        crtLoc = texAxes{k};
        
        % generate a texture all the way at the end of the texture axis
        generator = PatchAxisGenerator(plane, texAxes{k}, patchSize);
        generator.locations = 1;
        generator.next;
        crtPatch = generator.samples;
        
        % draw the generated patch
        fig = figure;
        fig.Units = 'pixels';
        fig.Position = [fig.Position(1:2) patchSize patchSize];
        
        ax = axes;
        ax.Position = [0 0 1 1];
        imagesc(crtPatch, [0 1]);
        colormap('gray');
        
        axis equal;
        axis off;
        
        % save the frame data as a PNG
        drawnow;
        frame = getframe(fig);
        imwrite(frame.cdata, fullfile('figs', 'draft', ...
            ['ternaryPatch_' plane '_' locLabels{k} '.png']));
    end
end

%% Draw simplex and guiding circles

% where to draw inset bar plots
texAxes = {[0 1 0], [0 0 1], [1 0 0]};
% color of bars
barColor = [64, 90, 128]/255;
% size of bar plot
patchDim = 0.5;
% shift of bar plot compared to axis coordinates
patchShift = 0.07;

% build the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2 2];

hold on;
axis equal;

% draw the guide simplex and circles
drawTernaryTriangle('edgelabels', 'none');

% draw
for k = 1:length(texAxes)
    crtLoc = texAxes{k};
    
    % project 3d axes coordinates to 2d alloy plot coordinates
    crtT = ternary3to2(crtLoc);
    
    if k == 3
        crtT(1) = crtT(1) - 4*patchShift;
    elseif k == 2
        crtT(2) = crtT(2) - patchDim;
        crtT(1) = crtT(1) - 4*patchShift;
    else
        crtT(2) = crtT(2) - patchDim;
    end
    
    % draw the bar plots
    % the weird ifs just position the plots in more-or-less reasonable places
    if k ~= 3
        insetbar([crtT + [0 1.15*patchDim] patchDim patchDim/2], crtLoc + 0.05, ...
            'color', barColor);
    else
        insetbar([crtT + [0 -0.65*patchDim] patchDim patchDim/2], crtLoc + 0.05, ...
            'color', barColor);
    end
end

% draw a bar plot in the center
insetbar([0 patchShift patchDim patchDim/2], [1/3 1/3 1/3] + 0.05, ...
    'color', barColor);

% update the looks of the axes
beautifygraph('ticks', 'off', 'ticklabels', false, ...
    'titlesize', 10, 'titleweight', 'normal', 'noaxes', true, ...
    'fontscale', 0.9);
preparegraph;

safePrint(fullfile('figs', 'draft', 'exampleTernaryPlane'));

%% Choose an image region to work with

% use the first image from the Penn database
images = parseImageNameFile('PennNoSkyIndex.txt', 'NaturalImages');
img0 = loadLUMImage(images{1});

% convert to log, as we do during processing
img0log = logTransform(img0);
% contrast adaptation doesn't affect the processing, but improves the
% visualization
img = contrastAdapt(img0log(3*128:5*128-1, 4*128:6*128-1), 'meanFct', @median);

% show the image region we'll use
figure;
imagesc(img);
axis equal;
colormap('gray');

%% Load filter

filter0 = open(fullfile('filters', 'filter4x32.mat'));
filter = filter0.filter;

%% Draw the preprocessing pipeline

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
% we'll use this function to brighten the images, otherwise they end up too
% dark when displayed
imgBrighten = @(m) 0.5*m + 0.4;

% select whether to draw patches with border or without
dimage = @image;
% dimage = @bimage;

% original image
% set up the sizes of the edges to use around the image, and the patch sizes
topEdge = 24;
leftEdge = 24;
patchDist = 24;
panelDist = (fig.Position(3)*2 - 256*5 - 2*leftEdge - patchDist/2) / 4;
h1 = dimage(leftEdge, topEdge, imgBrighten(img), 'CDataMapping', 'scaled');
colormap('gray');

% down-sampled (block-averaged) image
imgSmall = blockAverage(img, 4);
leftSmall = leftEdge + size(img, 2) + panelDist;
h2 = dimage([leftSmall leftSmall + size(img, 2)-1], ...
    [topEdge topEdge + size(img, 1)-1], imgBrighten(imgSmall), 'CDataMapping', 'scaled');

% draw patchified image
leftPatches = leftSmall + size(img, 2) + panelDist;
for i = 0:1
    % figure out position of top edge
    if i == 0
        topPatch = topEdge - patchDist/2;
    else
        topPatch = topEdge + size(img, 1)/2 + patchDist/2;
    end
    for j = 0:1
        % figure out position of left edge
        if j == 0
            leftPatch = leftPatches - patchDist / 2;
        else
            leftPatch = leftPatches + size(img, 2)/2 + patchDist/2;
        end
        % draw
        crtPatch = imgSmall(1+i*32:(i+1)*32, 1+j*32:(j+1)*32);
        dimage([leftPatch leftPatch + size(img, 2)/2-1], ...
            [topPatch topPatch + size(img, 1)/2-1], imgBrighten(crtPatch), ...
            'CDataMapping', 'scaled');
    end
end

% draw whitened patches; make sure to re-contrast adapt
imgSmallPatchified = patchify(imgSmall, size(filter));
% make sure the range for imgFiltered is similar to img
imgFiltered = contrastAdapt(filterImage(imgSmallPatchified, filter), ...
    'meanFct', @median);
leftWhitened = leftPatches + size(img, 2) + panelDist;
for i = 0:1
    if i == 0
        topPatch = topEdge - patchDist/2;
    else
        topPatch = topEdge + size(img, 1)/2 + patchDist/2;
    end
    for j = 0:1
        if j == 0
            leftPatch = leftWhitened - patchDist / 2;
        else
            leftPatch = leftWhitened + size(img, 2)/2 + patchDist/2;
        end
        crtPatch = imgFiltered(:, :, 2*j + i + 1);
        dimage([leftPatch leftPatch + size(img, 2)/2-1], ...
            [topPatch topPatch + size(img, 1)/2-1], imgBrighten(crtPatch), ...
            'CDataMapping', 'scaled');
    end
end

% draw ternarized patches
imgTernarized = quantizeImage(equalizeImage(imgFiltered), 3);
leftTernarized = leftWhitened + size(img, 2) + panelDist;
for i = 0:1
    % figure out position of top edge
    if i == 0
        topPatch = topEdge - patchDist/2;
    else
        topPatch = topEdge + size(img, 1)/2 + patchDist/2;
    end
    for j = 0:1
        % figure out position of left edge
        if j == 0
            leftPatch = leftTernarized - patchDist / 2;
        else
            leftPatch = leftTernarized + size(img, 2)/2 + patchDist/2;
        end
        % draw
        crtPatch = imgTernarized(:, :, 2*j + i + 1);
        dimage([leftPatch leftPatch + size(img, 2)/2-1], ...
            [topPatch topPatch + size(img, 1)/2-1], crtPatch, ...
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
fname = fullfile('figs', 'draft', 'ternaryNIPipelineDescription.png');
figFrame = getframe(fig);
pixels = figFrame.cdata;
imwrite(pixels, fname, 'png');

%% Store a single patch

crtPatch = imgTernarized(1:32, 1:32);
% imagesc(crt_patch, [0 1]);
bimage(0, 0, crtPatch, 'CDataMapping', 'scaled');
set(gca, 'clim', [0 1]);
colormap('gray');
axis equal;
axis off;

preparegraph;

safePrint(fullfile('figs', 'draft', 'exampleTernaryPatch'));

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
[~, crtHisto] = analyzeTexture(crtPatch, 3);

% draw the histogram
fig = figure;
fig.Units = 'inches';
fig.Position = [2 2 8 4];

ax = axes;
ax.OuterPosition = [0 0.2 1 0.8];

crtColor = [0.5 0.7 1]/2;
bar(crtHisto, 0.6, 'edgecolor', 'none', 'facecolor', crtColor);
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

gliderSize = 2;
topEdge = 10;
for i = 1:81
    row = mod(i-1, 3);
    bimage(i - gliderSize/2, topEdge + row*gliderSize*1.5, gliders{i}, ...
        'CDataMapping', 'scaled', 'gridwidth', 1);
end
colormap('gray');

xlim([0 82]);
% ylim([0 ax2.Position(4)/ax2.Position(3)*82]);

ax2.CLim = [0 2];

axis equal;
axis off;

preparegraph;

safePrint(fullfile('figs', 'draft', 'exampleTernaryHistogram'));

%% Show ternarization CDF

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.1 0.7];

% get a CDF for the patch
crtPatchFiltered = imgFiltered(1:32, 1:32);
[patchHist, patchBins] = hist(crtPatchFiltered(:), 32);
patchCdf = cumsum(patchHist) / sum(patchHist);

% find the positions of the cutoffs
cutoffsY = [1/3, 2/3];
cutoffsX = zeros(1, 2);
cutoffsIdx = zeros(1, 2);
for i = 1:length(cutoffsY)
    crtIdx = find(patchCdf > cutoffsY(i), 1);
    crtFrac = (cutoffsY(i) - patchCdf(crtIdx-1)) / (patchCdf(crtIdx) - patchCdf(crtIdx-1));
    cutoffsX(i) = (1-crtFrac)*patchBins(crtIdx-1) + crtFrac*patchBins(crtIdx);
    cutoffsIdx(i) = crtIdx - 1;
end

% fill([patch_bins(1) cutoffs_x(1) cutoffs_x(1) patch_bins(cutoffs_idx(1):-1:1)], ...
%      [0 0 cutoffs_y(1) patch_cdf(cutoffs_idx(1):-1:1)], 'k');
% fill([cutoffs_x(1) cutoffs_x(2) cutoffs_x(2) patch_bins(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_x(1)], ...
%      [0 0 cutoffs_y(2) patch_cdf(cutoffs_idx(2):-1:cutoffs_idx(1)+1) cutoffs_y(1)], ...
%      [0.5 0.5 0.5]);

% draw the CDF regions mapping to each gray level
hold on;

fill([patchBins(1) patchBins(1) cutoffsX(1) patchBins(cutoffsIdx(1):-1:1)], ...
     [0 cutoffsY(1) cutoffsY(1) patchCdf(cutoffsIdx(1):-1:1)], 'k', 'edgecolor', 'none');
fill([patchBins(1) patchBins(1) cutoffsX(2) patchBins(cutoffsIdx(2):-1:cutoffsIdx(1)+1) cutoffsX(1)], ...
     [cutoffsY(1) cutoffsY(2) cutoffsY(2) patchCdf(cutoffsIdx(2):-1:cutoffsIdx(1)+1) cutoffsY(1)], [0.5 0.5 0.5], ...
     'edgecolor', 'none');

% draw guide lines for the cutoffs
for i = 1:length(cutoffsX)
    plot([cutoffsX(i) patchBins(end)], repmat(cutoffsY(i), 1, 2), 'k:');
end

% draw the CDF
plot(patchBins, patchCdf, 'color', [0    0.4470    0.7410]);

% fix the axes (ticks, labels, etc.)
xlim(patchBins([1 end]));
tl = get(gca, 'ticklength');
set(gca, 'xtick', [], 'ytick', [0 1/3 2/3 1], 'yticklabel', {'0', '1/3', '2/3', '1'}, ...
    'yminortick', 'off', 'ticklength', [0 tl(2)]);

xlabel('luminance');
ylabel('cdf');

% beautify and save
beautifygraph('fontscale', 0.5, 'ticksize', 12);
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternarizationCdf'));

%% Show ternarization PDF

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.1 0.7];

hold on;
% draw the histogram (different colors for black, gray, white levels)
bar(patchBins(1:cutoffsIdx(1)), patchHist(1:cutoffsIdx(1)), 'k', 'linewidth', 0.2);
bar(patchBins(cutoffsIdx(1)+1:cutoffsIdx(2)), patchHist(cutoffsIdx(1)+1:cutoffsIdx(2)), ...
    'facecolor', [0.5 0.5 0.5], 'linewidth', 0.2);
bar(patchBins(cutoffsIdx(2)+1:end), patchHist(cutoffsIdx(2)+1:end), 'w', 'linewidth', 0.2);

% fix the axes (ticks, labels, etc.)
xlabel('luminance');
ylabel('count');

tl = get(gca, 'ticklength');
set(gca, 'xtick', [], 'ytick', max(patchHist)*[0 1/3 2/3 1], 'yticklabel', {'0', '1/3', '2/3', '1'}, ...
    'yminortick', 'off', 'ticklength', [0 tl(2)]);
xlim(patchBins([1 end]));
ylim([0 max(patchHist)]);

% beautify and save
beautifygraph('fontscale', 0.5, 'ticksize', 12);
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternarizationPdf'));