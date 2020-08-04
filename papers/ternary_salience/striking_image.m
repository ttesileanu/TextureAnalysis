% Making a striking image for eLife.

%% Load some data and predictions

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));
load(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32_square.mat'));

%% Choose the data

plane_choice = 'AB_1_2';

pp_mask = strcmp(ternaryAvg.groups, plane_choice);
pp_choice.thresholds = ternaryAvg.thresholds(pp_mask);
pp_choice.directions = ternaryAvg.directions(pp_mask);

ni_mask = strcmp(predictions.groups, plane_choice);
ni_choice.thresholds = predictions.thresholds(ni_mask);
ni_choice.directions = predictions.directions(ni_mask);

%% The hard part: make the ring of texture (or a disk, in this case)

% keep this reproducible
rng default;

% at least 1800 x 900 in size
total_size = 1200;
texture_image = zeros(total_size, total_size);

disk_width = 0.2;

crt_x_coords = repmat(linspace(-1 - disk_width, 1 + disk_width, total_size), ...
    total_size, 1);
crt_y_coords = repmat(linspace(-1 - disk_width, 1 + disk_width, total_size)', ...
    1, total_size);

crt_r_coords = sqrt(crt_x_coords .^2 + crt_y_coords .^ 2);
crt_theta_coords = atan2(crt_y_coords, crt_x_coords);
crt_theta_coords(crt_theta_coords < 0) = crt_theta_coords(crt_theta_coords < 0) + 2 * pi;

n_patches = 36;
patch_size = [32, 32];
scaled_size = [160, 160];
angle_step = 2 * pi / n_patches;
patch_angles = (0:n_patches - 1) * angle_step;
angle_space = angle_step / 10;
for i = 1:n_patches
    crt_angle = patch_angles(i);
    crt_x = cos(crt_angle);
    crt_y = sin(crt_angle);
    
    crt_loc = ternary2to3([crt_x, crt_y]);
    generator = PatchAxisGenerator(plane_choice, crt_loc, patch_size);
    generator.nLocations = 2;
    
    % this is a little silly, but it allows us to use PatchAxisGenerator's
    % code for figuring out where the maximum saturation is
    % skipping the unbiased sample
    generator.next;
    generator.next;
    crt_patch = generator.samples(1);
    
    % adapt colors
    crt_patch = 2 * crt_patch + 1;
    
    % scale up
    crt_patch_big = kron(crt_patch, ones(scaled_size ./ patch_size));
    
    % place the patch on the image
    angle_step_red = angle_step - angle_space;
    if crt_angle >= angle_step_red / 2 && crt_angle + angle_step_red / 2 <= 2 * pi
        crt_mask = (crt_theta_coords >= crt_angle - angle_step_red / 2) & ...
            (crt_theta_coords < crt_angle + angle_step_red / 2);
    elseif crt_angle < angle_step_red / 2
        crt_mask1 = crt_theta_coords < crt_angle + angle_step_red / 2;
        crt_mask2 = crt_theta_coords >= 2 * pi + crt_angle - angle_step_red / 2;
        crt_mask = crt_mask1 | crt_mask2;
    else  % crt_angle + angle_step / 2 > 2 * pi
        crt_mask1 = crt_theta_coords >= crt_angle - angle_step_red / 2;
        crt_mask2 = crt_theta_coords < crt_angle + angle_step_red / 2 - 2 * pi;
        crt_mask = crt_mask1 | crt_mask2;
    end
    
    crt_mask = crt_mask & (crt_r_coords >= 1) & (crt_r_coords <= 1 + disk_width);
    
    [mask_i, mask_j] = ind2sub(size(crt_mask), find(crt_mask));
    n_i = max(mask_i) - min(mask_i) + 1;
    n_j = max(mask_j) - min(mask_j) + 1;
    crt_subpatch_big = crt_patch_big(1 : n_i, 1 : n_j);
    
    texture_image(crt_mask) = 0;
    texture_image(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) = ...
        texture_image(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) + ...
        crt_mask(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) .* ...
        crt_subpatch_big;
end

% place unbiased random texture in the center
center_radius = 0.12;
center_patch = 1 + randi([0, 2], patch_size);
center_patch_big = kron(center_patch, ones(scaled_size ./ patch_size));

center_mask = crt_r_coords < center_radius;
texture_image(center_mask) = 0;

[mask_i, mask_j] = ind2sub(size(center_mask), find(center_mask));
n_i = max(mask_i) - min(mask_i) + 1;
n_j = max(mask_j) - min(mask_j) + 1;
crt_subpatch_big = center_patch_big(1 : n_i, 1 : n_j);
texture_image(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) = ...
        texture_image(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) + ...
        center_mask(min(mask_i) : max(mask_i), min(mask_j) : max(mask_j)) .* ...
        crt_subpatch_big;

%% Make a basic plot

% start plotting
[~, color_dict] = get_palette;

fig = figure;
fig.Units = 'inches';
fig.Position = [1, 2, 6, 4];
fig.InvertHardcopy = 'off';

ax = axes;
ax.Position = [0, 0, 1, 1];

hold on;

image([-1 - disk_width, 1 + disk_width], [-1 - disk_width, 1 + disk_width], ...
    1 + texture_image);
% bkg_color = color_dict('dark blue');
bkg_color = [1 1 1] / 255;
special_cmap = [bkg_color ; 0 0 0 ; 0.5 0.5 0.5 ; 0.99 0.99 0.99];
colormap(special_cmap);

% easy part: show the thresholds
plot_options = {'ellipse', true, 'ellipseOptions', {'nPoints', 128}};
plotTernaryThresholds(pp_choice, 'marker', 'x', 'color', color_dict('red'), ...
    plot_options{:});
plotTernaryThresholds(ni_choice, 'marker', '+', 'color', color_dict('orange'), ...
    plot_options{:});

ax.Color = bkg_color;

t = linspace(0, 2 * pi, 128);
x = cos(t);
y = sin(t);
plot(x, y, 'color', [0.5, 0.5, 0.5]);

% arrange axes and prepare for saving
axis equal;

xticks([]);
yticks([]);

xlim([-1.4, 2.7]);

beautifygraph;
preparegraph('edge', 0);

colormap(special_cmap);

% print('-dpng', '-r300', fullfile('figs', 'draft', 'striking_image.png'));
print('-dpdf', fullfile('figs', 'draft', 'striking_image.pdf'));

%% Store some natural-scene frames

img_names = {'DSC_0002_LUM.mat', 'DSC_0003_LUM.mat', 'DSC_0004_LUM.mat'};
path = fullfile('NaturalImages', 'PennNoSky', 'cd01A');

for i = 1:length(img_names)
    crt_path = fullfile(path, img_names{i});
    crt_img = loadLUMImage(crt_path);
    
    [ny, nx] = size(crt_img);
    
    crt_img = crt_img(floor(ny / 3) : end, ...
        floor(nx / 4) : floor(11 * nx / 12));
   
    crt_img = equalizeImage(crt_img);
    
    fig = figure;
    fig.Units = 'pixels';
    fig.Position = [fig.Position(1:2) fliplr(size(crt_img))];
    
    ax = axes;
    ax.Position =[0, 0, 1, 1];
    imagesc(crt_img);
    colormap('gray');
    axis equal;
    
    [~, crt_name_no_ext] = fileparts(img_names{i});
    
    out_fname = fullfile('figs', 'draft', [crt_name_no_ext '.png']);
    fig_frame = getframe(fig);
    crt_pixels = fig_frame.cdata;
    imwrite(crt_pixels, out_fname, 'png');
end
