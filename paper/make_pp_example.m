% make an example psychophysics trial

%% draw ternary strips

rng(34546);

patch_size = 64;
strip_size = 16;
strip_distance = 8;

bkg_patch = randi([0 2], patch_size)/2;

generator = PatchAxisGenerator('AD_1_2', [1.0 0.0 0.0], [patch_size strip_size]);
generator.locations = [0.2 0.4 0.7];
i = 1;
while generator.next
    strip = generator.samples;
    
    overlay_patch = bkg_patch;
    overlay_patch(:, strip_distance:strip_distance+strip_size-1) = strip;
    
    imagesc(overlay_patch);
    axis equal;
    colormap('gray');
    xlim([0.5 patch_size+0.5]);
    ylim([0.5 patch_size+0.5]);
    box on;
    set(gca, 'xtick', [], 'ytick', []);
    
    safe_print(fullfile('figs', 'draft', ['ternary_pp_example' int2str(i)]), ...
        'type', 'png');
    
    i = i + 1;
end

imagesc(bkg_patch);
axis equal;
colormap('gray');
xlim([0.5 patch_size+0.5]);
ylim([0.5 patch_size+0.5]);
box on;
set(gca, 'xtick', [], 'ytick', []);

safe_print(fullfile('figs', 'draft', 'ternary_pp_example_bkg'), 'type', 'png');

%% draw example psychophysics curve

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 1.5 1];

xvals = linspace(0, 1);
% yvals = 1./(1 + exp(-10*(xvals - 0.5)));
yvals = 1 - 2.^(-(xvals/0.5).^2.5);
color = [0.8 0.3 0.3];
plot(xvals, yvals, 'color', color);
hold on;
plot(xvals(end/2), yvals(end/2), '.', 'color', color);

set(gca, 'ytick', [0 1], 'yticklabel', {'chance', 'perfect'}, 'xtick', [0, 1], ...
    'xticklabel', {'1/3', '1'});
% xlabel('coordinate');

box off;

preparegraph;

safe_print(fullfile('figs', 'draft', 'example_pp_curve'));

%% draw ternary patches

rng(34546);

patch_size = 64;

generator = PatchAxisGenerator('AD_1_2', [1.0 0.0 0.0], [patch_size patch_size]);
generator.locations = [0.0 0.5 1.0];
i = 1;
while generator.next
    patch = generator.samples;
    
    fig = figure;
    fig.Units = 'pixels';
    fig.Position = [fig.Position(1:2) 2*patch_size 2*patch_size];
    
    ax = axes;
    ax.Position = [0 0 1 1];
    imagesc(patch, [0 1]);
    colormap('gray');
    
    axis equal;
    axis off;
    
    drawnow;
    frame = getframe(fig);
    imwrite(frame.cdata, fullfile('figs', 'draft', ...
        ['ternary_pp_full_patch' int2str(i) '.png']));
    
    i = i + 1;
end

%% show measurement rays, single plane

ternary_avg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 2 2];

drawTernaryTriangle('edgelabels', 'digit');
mask = strcmp(ternary_avg.groups, 'AD_1_2');
directions = ternary_avg.directions(mask);

axis equal;
axis off;

hold on;

color = [0.8 0.3 0.3];
for i = 1:length(directions)
    crt_dir = directions{i};
    % locate the position on the simplex
    min_coord = min(crt_dir(:));
    crt_dir = crt_dir / (1 - 3*min_coord);                        
    crt_dir2 = ternary3to2(crt_dir);
    plot([0 crt_dir2(1)], [0 crt_dir2(2)], 'color', color);
end

preparegraph;

safe_print(fullfile('figs', 'draft', 'single_plane_pp_directions'));

%% show measurement rays, mixed plane

ternary_avg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 2 2];

drawTernaryMixedBackground('Direction 1', 'Direction 2');
mask = strcmp(ternary_avg.groups, 'AB_1_1[2];AC_1_1[0]');
directions = cell2mat(ternary_avg.directions(mask));
directions2 = ternary6tomix2(directions);

axis equal;
axis off;

hold on;

color = [0.8 0.3 0.3];
for i = 1:length(directions)
    crt_dir2 = directions2(i, :);
    plot([0 crt_dir2(1)], [0 crt_dir2(2)], 'color', color);
end

preparegraph;

safe_print(fullfile('figs', 'draft', 'mixed_plane_pp_directions'));