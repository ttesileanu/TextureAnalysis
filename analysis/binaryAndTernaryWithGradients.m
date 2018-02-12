% try to fit binary and ternary data using gradient measures

%% Select blockAF and patch sizes

% [N, R] pairs
% NR_values = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
%     [4, 32], [4, 48], [4, 64]};
NR_values = {[2, 32]};

N_values = cellfun(@(x) x(1), NR_values);
R_values = cellfun(@(x) x(2), NR_values);

%% Load multiscale filters

filters = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
    filterFilename = fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']);
    crtFilter = open(filterFilename);
    fields = fieldnames(crtFilter);
    isfilter = @(f) isnumeric(f) && ismatrix(f) && all(size(f) == [crtR crtR]);
    valid = cellfun(@(s) isfilter(crtFilter.(s)), fields);
    if sum(valid) == 0
        error(['Can''t find valid filter data in ' filterFilename '.']);
    elseif sum(valid) > 1
        error(['Don''t know which field to use from ' filterFilename '.']);
    end
    filterField = fields{valid};
    filters{i} = crtFilter.(filterField);
end

%% Generate the gradient stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_gradients = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on gradient stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_gradients{i} = walkImageSet(@walkerGradients, images, 'NaturalImages', ...
        'filter', filters{i}, 'equalize', 'contrast', ...
        'args', {'analysisPatchSize', crtR});
end

%% Save the continuous stats

res = res_gradients;
save(fullfile('save', 'natural_nosky_gradients_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Choose a value for N and R

gradient_stats = res_gradients{1};

%% Analyze the distribution using histograms

nev = size(gradient_stats.ev, 2);
nx = ceil(sqrt(nev));
ny = ceil(nev/nx);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 3*nx 2*ny];
for i = 1:nev
    subplot(ny, nx, i);
    hist(gradient_stats.ev(:, i), 100);
    crtLabel = gradient_stats.statsDesc{i};
    crtLabel(crtLabel == '_') = ' ';
    xlabel(crtLabel);
    beautifygraph;
end

preparegraph;

%% Analyze the distribution using scatter plots

plot_pairs = {[1, 2], [1, 3], [2, 4], [1, 5], [1, 6], [5, 7], [5, 8], [6, 7], [1, 8]};
npairs = length(plot_pairs);
nx = ceil(sqrt(npairs));
ny = ceil(nev/nx);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 3*nx 2.5*ny];
for i = 1:npairs
    subplot(ny, nx, i);
    
    i1 = plot_pairs{i}(1);
    i2 = plot_pairs{i}(2);
    
    smartscatter(gradient_stats.ev(:, i1), gradient_stats.ev(:, i2));
    
    crtLabel1 = gradient_stats.statsDesc{i1};
    crtLabel1(crtLabel1 == '_') = ' ';
    crtLabel2 = gradient_stats.statsDesc{i2};
    crtLabel2(crtLabel2 == '_') = ' ';
    
    xlabel(crtLabel1);
    ylabel(crtLabel2);
    
    beautifygraph;
end

preparegraph;

%% Analyze the distribution using PCA

[all_pc_coeffs, all_pc_proj, all_pc_vars, ~, all_pc_explained] = pca(gradient_stats.ev);

labels = gradient_stats.statsDesc;
for i = 1:length(labels)
    labels{i}(labels{i} == '_') = ' ';
end

% show top EVs
nev = size(gradient_stats.ev, 2);
n_pcs = nev;
fig = figure;
fig.Units = 'inches';
fig.Position = [0 1 2.5*n_pcs 3];
for i = 1:n_pcs
    crt_pc = all_pc_coeffs(:, i);
    subplot(1, n_pcs, i);
    
    bar(crt_pc);
    
    xlim([0.5, 0.5 + nev]);
    ylim([-1 1]);
    
    title(['PC' int2str(i)]);
    
    set(gca, 'xtick', 1:nev, 'xticklabel', labels, 'xticklabelrotation', 90);
end

% show distribution of natural images in space of top PCs
figure;
smartscatter(all_pc_proj(:, 1), all_pc_proj(:, 2));
xlabel('PC1');
ylabel('PC2');

%% Converting from binary patches to gradient stats

% how many locations along each axis
n_locs_binary = 8;

% one axis in each of the coordinate directions, but each must appear twice
% (for + and - directions)!
binary_map_groups0 = {'A_1', 'AC_1_1', 'AB_1_1', 'AD_1_1', 'BC_1_1', ...
    'ABC_1_1_1', 'BCD_1_1_1', 'ABD_1_1_1', 'ACD_1_1_1', 'ABCD_1_1_1_1'};
binary_map_groups = flatten(repmat(binary_map_groups0, 2, 1));
binary_map_axes = flatten([repmat({[1 0]}, 1, length(binary_map_groups0)) ; ...
    repmat({[0 1]}, 1, length(binary_map_groups0))]);

binary_map_locs = cell(size(binary_map_groups));

% number of patches to generate at each location
n_samples_binary = 16;
binary_map_patches = cell(size(binary_map_groups));
binary_map_stats = cell(size(binary_map_groups));

% size of patch
patch_size = 64;

nev = size(gradient_stats.ev, 2);
progress = TextProgress('generating binary patches', 'prespace', 32, 'length', 10);
for i = 1:length(binary_map_groups)
    % generate patches in this direction
    generator = PatchAxisGenerator(binary_map_groups{i}, binary_map_axes{i}, ...
        patch_size);
    generator.nLocations = n_locs_binary;
    
    crt_patches = cell(1, n_locs_binary);
    j = 1;
    while generator.next
        crt_patches{j} = generator.samples(n_samples_binary);
        j = j + 1;
    end
    % store patches
    binary_map_patches{i} = crt_patches;
    
    % process patches with Chubb nonlinearity, then calculate continuous stats
    crt_evs = cell(1, n_locs_binary);
    for j = 1:n_locs_binary
        loc_ev = zeros(n_samples_binary, nev);
        for p = 1:n_samples_binary
            crt_patch = crt_patches{j}(:, :, p);
            crt_grad_res = analyzePatchGradients(crt_patch, []);
            loc_ev(p, :) = crt_grad_res.ev;
        end
        crt_evs{j} = loc_ev;
    end
    % store stats
    binary_map_stats{i} = crt_evs;
    % store locations
    binary_map_locs{i} = generator.getLocations;
    
    progress.update(100*i/length(binary_map_groups));
end
progress.done('done');

%% Load binary psychophysics data

binary_data0 = open(fullfile('data', 'PPandNIstats.mat'));
binary_thresholds = 1./sqrt(1000*diag(binary_data0.dataPP.allS.subjAvg.covM));

%% Load ternary psychophysics data

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

%% Converting from ternary patches to gradient stats

% how many locations along each axis
n_locs_ternary = 8;

% one axis in each of the coordinate directions
ternary_map_groups = ternary_groups;
ternary_map_locs = cell(size(ternary_map_groups));
ternary_map_axes = ternary_axes;

% number of patches to generate at each location
n_samples_ternary = 16;
ternary_map_patches = cell(size(ternary_map_groups));
ternary_map_stats = cell(size(ternary_map_groups));

% size of patch
patch_size = 64;

nev = size(gradient_stats.ev, 2);
progress = TextProgress('generating ternary patches', 'prespace', 28, 'length', 20);
for i = 1:length(ternary_map_groups)
    % generate patches in this direction
    generator = PatchAxisGenerator(ternary_map_groups{i}, ternary_map_axes{i}, ...
        patch_size);
    generator.nLocations = n_locs_ternary;
    
    crt_patches = cell(1, n_locs_ternary);
    j = 1;
    while generator.next
        crt_patches{j} = generator.samples(n_samples_ternary);
        j = j + 1;
    end
    % store patches
    ternary_map_patches{i} = crt_patches;
    
    % process patches with Chubb nonlinearity, then calculate continuous stats
    crt_evs = cell(1, n_locs_ternary);
    for j = 1:n_locs_ternary
        loc_ev = zeros(n_samples_ternary, nev);
        for p = 1:n_samples_ternary
            crt_patch = crt_patches{j}(:, :, p);
            crt_grad_res = analyzePatchGradients(crt_patch, []);
            loc_ev(p, :) = crt_grad_res.ev;
        end
        crt_evs{j} = loc_ev;
    end
    % store stats
    ternary_map_stats{i} = crt_evs;
    % store locations
    ternary_map_locs{i} = generator.getLocations;
    
    progress.update(100*i/length(ternary_map_groups));
end
progress.done('done');

%% Find location of unbiased random patches

test_patch_size = 64;
test_nsamples = 1024;
nev = size(gradient_stats.ev, 2);
test_random_stats = zeros(test_nsamples, nev);

for k = 1:test_nsamples
    crt_patch = rand(test_patch_size);
    crt_grad_res = analyzePatchGradients(crt_patch, []);
    test_random_stats(k, :) = crt_grad_res.ev;
end

origin = mean(test_random_stats, 1);

%% Plot coordinate axes in PC1/2 space

% calculate PCA for natural image patches
all_pc_coeffs = pca(gradient_stats.ev);

ni_all_ev_proj_nonmc = gradient_stats.ev*all_pc_coeffs;

figure;
hold on;
smartscatter(ni_all_ev_proj_nonmc(:, 1), ni_all_ev_proj_nonmc(:, 2), 'alpha', 0.1);

% bin_axes_to_show = 1:length(binary_map_groups);
% focus on only second-order axes
bin_axes_to_show = find(cellfun(@(s) length(s) == 6, binary_map_groups));
binary_interpolated = mapInterpolate(binary_map_stats, binary_map_locs, 2);
trajectory_locs = linspace(0, 1, 10);
for i0 = 1:length(bin_axes_to_show)
    i = bin_axes_to_show(i0);
    crt_trajectory = binary_interpolated{i}.function(trajectory_locs);
    crt_trajectory_proj = crt_trajectory*all_pc_coeffs;
    plot(crt_trajectory_proj(:, 1), crt_trajectory_proj(:, 2), 'k');
    text(crt_trajectory_proj(end, 1), crt_trajectory_proj(end, 2), binary_map_groups{i});
end

ternary_interpolated = mapInterpolate(ternary_map_stats, ternary_map_locs, 2);
% focus on only second-order axes
tern_axes_to_show = find(cellfun(@(s) length(s) == 6, ternary_map_groups));
%tern_axes_to_show = 1:length(ternary_map_groups);
for i0 = 1:length(tern_axes_to_show)
    i = tern_axes_to_show(i0);
    crt_trajectory = ternary_interpolated{i}.function(trajectory_locs);
    crt_trajectory_proj = crt_trajectory*all_pc_coeffs;
    plot(crt_trajectory_proj(:, 1), crt_trajectory_proj(:, 2), 'b');
    text(crt_trajectory_proj(end, 1), crt_trajectory_proj(end, 2), ternary_map_groups{i});
end

origin_projected = origin(:)'*all_pc_coeffs;
plot(origin_projected(1), origin_projected(2), 'kx', 'linewidth', 2);

xlabel('PC1');
ylabel('PC2');

beautifygraph;
preparegraph;

%% Plot the measured binary and ternary thresholds in PC1/2 space

% find locations of binary thresholds
binary_interpolated = mapInterpolate(binary_map_stats, binary_map_locs, 2);
binary_threshold_locs = arrayfun(@(i) ...
    binary_interpolated{i}.function(binary_thresholds(i)), ...
    1:length(binary_thresholds), 'uniform', false);

% find locations of ternary thresholds
ternary_interpolated = mapInterpolate(ternary_map_stats, ternary_map_locs, 2);
ternary_threshold_locs = arrayfun(@(i) ...
    ternary_interpolated{i}.function(ternary_thresholds(i)), ...
    1:length(ternary_thresholds), 'uniform', false);

% calculate PCA for natural image patches
all_pc_coeffs = pca(gradient_stats.ev);

% project binary and ternary thresholds
binary_threshold_locs_projected = cellfun(@(v) ...
    v(:)'*all_pc_coeffs, binary_threshold_locs, 'uniform', false);
ternary_threshold_locs_projected = cellfun(@(v) ...
    v(:)'*all_pc_coeffs, ternary_threshold_locs, 'uniform', false);

% show distribution of natural images in space of top PCs
figure;
ni_all_ev_proj_nonmc = gradient_stats.ev*all_pc_coeffs;
smartscatter(ni_all_ev_proj_nonmc(:, 1), ni_all_ev_proj_nonmc(:, 2));
% overlap the threshold positions
hold on;
binary_threshold_locs_projected_mat = cell2mat(binary_threshold_locs_projected');
h_bin = scatter(binary_threshold_locs_projected_mat(:, 1), ...
    binary_threshold_locs_projected_mat(:, 2), [], [0.8 0 0], '.');
ternary_threshold_locs_projected_mat = cell2mat(ternary_threshold_locs_projected');
h_tern = scatter(ternary_threshold_locs_projected_mat(:, 1), ...
    ternary_threshold_locs_projected_mat(:, 2), [], [0 0.7 0], '.');

origin_projected = origin(:)'*all_pc_coeffs;
plot(origin_projected(1), origin_projected(2), 'kx', 'linewidth', 2);

xlabel('PC1');
ylabel('PC2');

legend([h_bin h_tern], {'binary', 'ternary'});

beautifygraph;
preparegraph;

%% Find threshold predictions

ni_cov = cov(gradient_stats.ev);
binary_thresholds_rep = flatten(repmat(binary_thresholds(:)', 2, 1));

trafo_matrix = sqrtm(ni_cov);
[binary_predictions, binary_mapped, binary_optim_details, binary_intersect_details] = ...
    optimizePredictedThresholds(binary_thresholds_rep, binary_interpolated, ...
    trafo_matrix, [0 0.5], 'ignore_mask', strcmp(binary_map_groups, 'A_1'));

[ternary_predictions, ternary_mapped, ternary_optim_details, ternary_intersect_details] = ...
    optimizePredictedThresholds(ternary_thresholds, ternary_interpolated, ...
    trafo_matrix, [0 0.5], 'ignore_mask', strcmp(ternary_map_groups, 'A_1'));