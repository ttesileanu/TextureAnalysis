% study the effects of noise on texture results

%% Generate random samples, G=3 and continuous stats

rng(3443);

G = 3;
n_samples = 1000;
patch_size = 64;
% patches = arrayfun(@(p) randi([0, G], patch_size)/(G-1), 1:n_samples, 'uniform', false);
patches = arrayfun(@(p) rand(patch_size), 1:n_samples, 'uniform', false);
evsFiniteG = [];
evsCont = [];
progress = TextProgress;
for i = 1:n_samples
    [~, crt_ev] = processBlock(quantize(patches{i}, G), G);
    [~, crt_ev_cont] = processBlock(patches{i}, inf);
    if isempty(evsFiniteG)
        evsFiniteG = zeros(n_samples, length(crt_ev));
    end
    if isempty(evsCont)
        evsCont = zeros(n_samples, length(crt_ev_cont));
    end
    evsFiniteG(i, :) = crt_ev; %#ok<SAGROW>
    evsCont(i, :) = crt_ev_cont; %#ok<SAGROW>
    if mod(i, 10) == 0
        progress.update(100*i/n_samples);
    end
end
progress.done;

%% Shape of noise around origin, G=3 stats

% keep things reproducible
rng(32423);

% find all group names
mtc = processBlock('mtc', G);
groups = mtc.coord_groups;
directions = cell(1, G-1);
for i = 1:G-1
    crt_dir = num2str((1:G) == i, '%d');
    directions{i} = ['[' crt_dir ']'];
end

% generate labels for all the coordinates in ev
labels = cell(1, length(groups)*(G-1));
for i = 1:length(labels)
    crt_group_idx = floor((i-1) / (G-1)) + 1;
    crt_dir_idx = mod(i-1, G-1) + 1;
    labels{i} = [groups{crt_group_idx}.name directions{crt_dir_idx}];
end

% choose how many plots to make
n_coords = size(evsFiniteG, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFiniteG(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

ymax = 0;

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFiniteG(:, i)) + hist_bins;
    histogram(evsFiniteG(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
    
    yl = ylim;
    ymax = max(ymax, yl(2));
end

ax = plotter.get_axes();
for i = 1:length(ax)
    ylim(ax(i), [0 ymax]);
end

%% Shape of noise around origin, continuous stats

% keep things reproducible
rng(32423);

% find all group names
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% choose how many plots to make
n_coords = size(evsCont, 2);
n_plots = n_coords;

% choose a range for the histograms
group_stds = arrayfun(@(i) std(evsCont(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

ymax = 0;

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsCont(:, i)) + hist_bins;
    histogram(evsCont(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
    
    yl = ylim;
    ymax = max(ymax, yl(2));
end

ax = plotter.get_axes();
for i = 1:length(ax)
    ylim(ax(i), [0 ymax]);
end

%% Susceptibility to noise away from origin

rng(23498);

patch_size = 64;
% generator = PatchAxisGenerator('AB_1_1', [0 0 1], patch_size);
generator = PatchAxisGenerator('AB_1_1', [0.1140 0.9302 -0.0442], patch_size);
% generator = PatchAxisGenerator('AB_1_2', [0 0 1], patch_size);
generator.locations = 0.5;
generator.next;
base_patch = generator.samples;

G = 3;
n_samples = 1000;
noise_size = 0.3;
patchesFar = arrayfun(@(p) min(max(base_patch + noise_size*(2*rand(patch_size)-1), 0), 1), ...
    1:n_samples, 'uniform', false);
evsFarFiniteG = [];
evsFarCont = [];
progress = TextProgress;
for i = 1:n_samples
    [~, crt_ev] = processBlock(quantize(patchesFar{i}, G), G);
    [~, crt_ev_cont] = processBlock(patchesFar{i}, inf);
    if isempty(evsFarFiniteG)
        evsFarFiniteG = zeros(n_samples, length(crt_ev));
    end
    if isempty(evsFarCont)
        evsFarCont = zeros(n_samples, length(crt_ev_cont));
    end
    evsFarFiniteG(i, :) = crt_ev; %#ok<SAGROW>
    evsFarCont(i, :) = crt_ev_cont; %#ok<SAGROW>
    if mod(i, 10) == 0
        progress.update(100*i/n_samples);
    end
end
progress.done;

%% Shape of noise far from origin, G=3 stats

% keep things reproducible
rng(32423);

% find all group names
mtc = processBlock('mtc', G);
groups = mtc.coord_groups;
directions = cell(1, G-1);
for i = 1:G-1
    crt_dir = num2str((1:G) == i, '%d');
    directions{i} = ['[' crt_dir ']'];
end

% generate labels for all the coordinates in ev
labels = cell(1, length(groups)*(G-1));
for i = 1:length(labels)
    crt_group_idx = floor((i-1) / (G-1)) + 1;
    crt_dir_idx = mod(i-1, G-1) + 1;
    labels{i} = [groups{crt_group_idx}.name directions{crt_dir_idx}];
end

% choose how many plots to make
n_coords = size(evsFarFiniteG, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFarFiniteG(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

ymax = 0;

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFarFiniteG(:, i)) + hist_bins;
    histogram(evsFarFiniteG(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
    
    yl = ylim;
    ymax = max(ymax, yl(2));
end

ax = plotter.get_axes();
for i = 1:length(ax)
    ylim(ax(i), [0 ymax]);
end

%% Shape of noise far from origin, continuous stats

% keep things reproducible
rng(32423);

% find all group names
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% choose how many plots to make
n_coords = size(evsFarCont, 2);
n_plots = n_coords;

% set a range for the histograms
group_stds = arrayfun(@(i) std(evsFarCont(:, i)), 1:n_coords);
hist_range = 2*max(group_stds);
hist_bins = linspace(-hist_range, hist_range, 20);

% make the plots
plotter = MatrixPlotter(n_plots);

ymax = 0;

while plotter.next
    i = plotter.index;
    bin_edges = mean(evsFarCont(:, i)) + hist_bins;
    histogram(evsFarCont(:, i), bin_edges);
    
    xlim([bin_edges(1), bin_edges(end)]);
    
    beautifygraph('fontscale', 0.5);
    
    xlabel(labels{i}, 'fontsize', 10);
    
    yl = ylim;
    ymax = max(ymax, yl(2));
end


ax = plotter.get_axes();
for i = 1:length(ax)
    ylim(ax(i), [0 ymax]);
end

%% Load the psychophysics data

ternary_data_ext0 = open(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% select single planes, full data
ternary_all_edirs = fieldnames(ternary_data_ext0.ds_merged.avg.edirs);
ternary_edir_mask = true(size(ternary_all_edirs));
for i = 1:length(ternary_all_edirs)
    crt_edir = ternary_all_edirs{i};
    crt_edir_data = ternary_data_ext0.ds_merged.avg.edirs.(crt_edir);
    
    crt_cgroup_names = crt_edir_data.cgroup_names;
    if ~strcmp(crt_cgroup_names{1}, crt_cgroup_names{2})
        % this is not single plane; skip
        ternary_edir_mask(i) = false;
        continue;
    end
    
    if ismember([crt_edir 'full'], ternary_all_edirs)
        % keep only full
        ternary_edir_mask(i) = false;
        continue;
    end

%     % keep only fields not ending in 'full'
%     if length(crt_edir) > 4 && strcmp(crt_edir(end-3:end), 'full')
%         edir_mask(i) = false;
%         continue;
%     end
    
    if ismember(crt_edir, {'A1G', 'A1Ginter'})
        % keep only merge
        ternary_edir_mask(i) = false;
        continue;
    end
end
ternary_sub_edirs = ternary_all_edirs(ternary_edir_mask);

% convert data to something similar to the CSV file
ternary_groups = {};
ternary_axes = {};
ternary_thresholds = [];
ternary_stdev_ratios = [];
ternary_use_mc = false;
for i = 1:length(ternary_sub_edirs)
    crt_edir = ternary_sub_edirs{i};
    crt_edir_data = ternary_data_ext0.ds_merged.avg.edirs.(crt_edir);
    crt_group = crt_edir_data.cgroup_names{1};
    crt_uvecs = crt_edir_data.uvecs;
    crt_ndirs = size(crt_uvecs, 1);
    if ~ternary_use_mc
        % use combined (_all) data
        % XXX should use _comth instead of _all!
        crt_thresh = crt_edir_data.thresh_mags_all;
        % XXX are these actually stdev ratios?
        crt_thresh_hi = crt_edir_data.thresh_mags_ebhi_all;
        crt_thresh_lo = crt_edir_data.thresh_mags_eblo_all;
    else
        % use data from MC
        crt_thresh = crt_edir_data.thresh_mags_mc;
        crt_thresh_hi = crt_edir_data.thresh_mags_ebhi_mc;
        crt_thresh_lo = crt_edir_data.thresh_mags_eblo_mc;
    end
    
    crt_stdev_ratios = crt_thresh_hi./crt_thresh_lo;
    
    ternary_groups = [ternary_groups repmat({crt_group}, 1, crt_ndirs)]; %#ok<AGROW>
    
    % transform 2d coordinates to 3d probabilities
    crt_axes = 1/3 + crt_uvecs*[[2/3 ; -1/3] [-1/3 ; 2/3] [-1/3 ; -1/3]];
    crt_axes_cell = mat2cell(crt_axes, ones(crt_ndirs, 1), 3);
    
    ternary_axes = [ternary_axes crt_axes_cell(:)']; %#ok<AGROW>
    
    ternary_thresholds = [ternary_thresholds crt_thresh(:)']; %#ok<AGROW>
    ternary_stdev_ratios = [ternary_stdev_ratios crt_stdev_ratios(:)']; %#ok<AGROW>
end

%% Generate noise around all of the axes in one plane

% keep things reproducible
% rng(23498);
rng(84594);

patch_size = 64;
axis_location = 0.5;

G = 3;
n_samples = 1000;
noise_size = 0.3;
% noise_size = 0;

plane = 'AB_1_1';
plane_mask = strcmp(ternary_groups, plane);
n_axes = sum(plane_mask);
noises = cell(1, n_axes);

masked_axes = ternary_axes(plane_mask);

for i = 1:length(noises)
    crt_axis = masked_axes{i};
   
    generator = PatchAxisGenerator(plane, crt_axis, patch_size);
    generator.locations = axis_location*sqrt(2)/sqrt(3*norm(crt_axis)^2 - 1);
    generator.next;
    base_patch = generator.samples;

    crt_patches = arrayfun(@(p) min(max(base_patch + noise_size*(2*rand(patch_size)-1), 0), 1), ...
        1:n_samples, 'uniform', false);
    crt_evs = [];
    progress = TextProgress([int2str(i) '/' int2str(length(noises))]);
    for j = 1:n_samples
        [~, crt_ev] = processBlock(quantize(crt_patches{j}, G), G);
        if isempty(crt_evs)
            crt_evs = zeros(n_samples, length(crt_ev));
        end
        crt_evs(j, :) = crt_ev; %#ok<SAGROW>
        if mod(j, 10) == 0
            progress.update(100*j/n_samples);
        end
    end
    progress.done;
    
    noises{i} = crt_evs;
end

%% Draw the noise

figure;

hold on;

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

% draw circles for orientation (radius 1 and 1/2)
plot(cos(angle_range), sin(angle_range), ':', 'color', [0.4 0.4 0.4]);
plot(0.5*cos(angle_range), 0.5*sin(angle_range), ':', 'color', [0.4 0.4 0.4]);

% draw the probability triangle
plot([-1/2 1 -1/2 -1/2], [-max_t2 0 max_t2 -max_t2], 'color', [0.5 0.7 1]);

% draw the main axes
plot([0 1.5], [0 0], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 -1/2*1.5], [0 1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);
plot([0 -1/2*1.5], [0 -1.5*max_t2], ':', 'color', [1 0.6 0.6], 'linewidth', 1);

% label the corners
text(1.07, -0.1, '[0,1,0]', 'fontsize', 12);
text(-1.1,  max_t2+0.05, '[0,0,1]', 'fontsize', 12);
text(-1.1, -max_t2-0.01, '[1,0,0]', 'fontsize', 12);

axis equal;

% get mapping from the 99 probabilities to texture groups
mtc = processBlock('mtc', 3);
unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);

idx1 = find(strcmp(unique_groups, plane), 1);

for i = 1:length(noises)
    crt_noise = noises{i};
    crt_noise_full = reshape(crt_noise, [], 2, 33);
    crt_noise_full(:, 3, :) = 1 - sum(crt_noise_full, 2);
    
    crt_locs = crt_noise_full(:, :, idx1);
    
    % this is a litle roundabout, but easier conceptually
    crt_t1 = (3*crt_locs(:, 2) - 1)/2;
    crt_t2 = (crt_locs(:, 3) - crt_locs(:, 1)) * max_t2;
    
    plot(crt_t1, crt_t2, 'k.', 'linewidth', 1)
end

beautifygraph;
preparegraph;