% make plots for an example psychophysics trial

%% draw ternary strips

% keep things reproducible
rng(34546);

% set the size of the background patch, size of the strip, and distance
% from strip to edge
patchSize = 64;
stripSize = 16;
stripDistance = 8;

% the background patch
bkgPatch = randi([0 2], patchSize)/2;

% generate a non-trivial texture, for the strip
generator = PatchAxisGenerator('AD_1_2', [1.0 0.0 0.0], [patchSize stripSize]);
% use three different strengths for the texture
generator.locations = [0.2 0.4 0.7];
i = 1;
while generator.next
    strip = generator.samples;
    
    % this is the patch with the strip on top of the background
    overlayPatch = bkgPatch;
    overlayPatch(:, stripDistance:stripDistance+stripSize-1) = strip;

    % draw and save
    imagesc(overlayPatch);
    axis equal;
    colormap('gray');
    xlim([0.5 patchSize+0.5]);
    ylim([0.5 patchSize+0.5]);
    box on;
    set(gca, 'xtick', [], 'ytick', []);
    
    safePrint(fullfile('figs', 'draft', ['ternaryPPExample' int2str(i)]), ...
        'type', 'png');
    
    i = i + 1;
end

% also save the background patch without any strip
imagesc(bkgPatch);
axis equal;
colormap('gray');
xlim([0.5 patchSize+0.5]);
ylim([0.5 patchSize+0.5]);
box on;
set(gca, 'xtick', [], 'ytick', []);

safePrint(fullfile('figs', 'draft', 'ternaryPPExampleBkg'), 'type', 'png');

%% draw example psychophysics curve

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.5 1];

% make an arbitrary Weibull curve
xVals = linspace(0, 1);
% yvals = 1./(1 + exp(-10*(xvals - 0.5)));
yVals = 1/4 + (3/4) * (1 - 2.^(-(xVals/0.5).^2.5));

% draw the curve
color = [0.8 0.3 0.3];
plot(xVals, yVals, 'color', color);
hold on;
plot(xVals(end/2), yVals(end/2), '.', 'color', color);

% set(gca, 'ytick', [0.25 1], 'yticklabel', {'chance', 'perfect'}, 'xtick', [0, 1], ...
%     'xticklabel', {'1/3', '1'});
set(gca, 'ytick', [0.25 1], 'yticklabel', {'0.25', '1.00'}, 'xtick', [0, 1], ...
    'xticklabel', {'1/3', '1'});
% xlabel('coordinate');

box off;

ylim([0, 1]);

preparegraph;

% save to file
safePrint(fullfile('figs', 'draft', 'examplePPCurve'));

%% draw ternary patches

% keep things reproducible
rng(34546);

patchSize = 64;

% generate patches at a few locations along the axis
generator = PatchAxisGenerator('AD_1_2', [1.0 0.0 0.0], [patchSize patchSize]);
generator.locations = [0.0 0.5 1.0];
i = 1;
while generator.next
    patch = generator.samples;
    
    % draw the patch
    fig = figure;
    fig.Units = 'pixels';
    fig.Position = [fig.Position(1:2) 2*patchSize 2*patchSize];
    
    ax = axes;
    ax.Position = [0 0 1 1];
    imagesc(patch, [0 1]);
    colormap('gray');
    
    axis equal;
    axis off;
    
    % save to file
    drawnow;
    frame = getframe(fig);
    imwrite(frame.cdata, fullfile('figs', 'draft', ...
        ['ternaryPPFullPatch' int2str(i) '.png']));
    
    i = i + 1;
end

%% show measurement rays, single plane

% only loading the data for the directions in which thresholds are measured
ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 2 2];

% draw the guiding simplex&circles
drawTernaryTriangle('edgelabels', 'digit');
% choose the directions in some particular plane, it doesn't matter which
mask = strcmp(ternaryAvg.groups, 'AD_1_2');
directions = ternaryAvg.directions(mask);

% remove the axes
axis equal;
axis off;

hold on;

% draw the rays along which we have measurements
color = [0.8 0.3 0.3];
for i = 1:length(directions)
    crtDir = directions{i};
    % locate the position on the simplex
    minCoord = min(crtDir(:));
    crtDir = crtDir / (1 - 3*minCoord);                        
    crtDir2 = ternary3to2(crtDir);
    plot([0 crtDir2(1)], [0 crtDir2(2)], 'color', color);
end

preparegraph;

% save to file
safePrint(fullfile('figs', 'draft', 'singlePlanePPDirections'));
