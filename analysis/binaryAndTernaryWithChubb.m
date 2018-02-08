% compare predictions from natural images analyzed using Chubb
% nonlinearities and the binary and ternary psychophysics

%% Load natural image data

ni_nonlinear = open(fullfile('save', 'natural_nosky_nonlinear_contrastadapt.mat'));
ni_choice = 1;

ni_all_ev0 = cell2mat(cellfun(@(s) s.ev, ni_nonlinear.res(ni_choice, :), 'uniform', false));

% clean the data
% (XXX some of the patches were probably too bright for contrast adapt? check! XXX)
cleaning_mask = (ni_all_ev0(:, 1) > -0.5);
ni_all_ev = ni_all_ev0(cleaning_mask, :);

%% Load Chubb nonlinearities and find origin of coordinate system

chubb_nonlin = struct;
chubb_nonlin.pvals = cell2mat(cellfun(@(r) r.nonlinearity, ni_nonlinear.res, 'uniform', false));

% Origin is where i.i.d. images with uniform histogram would map.
% For an n-point function, this should be (2*m - 1)^n, where m is the mean
% intensity value.
n_nonlin = size(chubb_nonlin.pvals, 2);
origin = zeros(1, 10*n_nonlin);
for i = 1:n_nonlin
    crt_base = (i-1)*10 + 1;
    crt_mean_scaled = 2*mean(chubb_nonlin.pvals(:, i)) - 1;
    origin(crt_base) = crt_mean_scaled;
    origin(crt_base+1:crt_base+4) = crt_mean_scaled^2;
    origin(crt_base+5:crt_base+8) = crt_mean_scaled^3;
    origin(crt_base+9) = crt_mean_scaled^4;
end

%% Check that the origin is correct by generating random samples

% there will be small discrepancies because of edge effects
test_patch_size = 64;
test_nsamples = 1024;
test_random_stats = zeros(test_nsamples, 10, n_nonlin);

for k = 1:test_nsamples
    crt_patch = rand(test_patch_size);
    for i = 1:n_nonlin
        crt_patch_nonlin = applyNonlinearity(crt_patch, chubb_nonlin.pvals(:, i));
        [~, crt_ev] = processBlock(crt_patch_nonlin, inf);
        test_random_stats(k, :, i) = crt_ev;
    end
end

test_random_stats_mean = flatten(mean(test_random_stats, 1));

%% Analyze the distribution using PCA

[all_pc_coeffs, all_pc_proj, all_pc_vars, ~, all_pc_explained] = pca(ni_all_ev);

% show top EVs
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
          '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

n_pcs = 10;
fig = figure;
fig.Units = 'inches';
fig.Position = [0 1 15 6];
crt = 1;
for j = 1:n_nonlin
    for i = 1:n_pcs
        crt_pc = all_pc_coeffs(:, i);
        subplot(n_nonlin, n_pcs, crt);
                
        bar(crt_pc((j-1)*10 + 1:j*10));
        
        crt = crt + 1;
        
        xlim([0.5, 10.5]);
        ylim([-1 1]);
        
        if j == 1
            title(['PC' int2str(i)]);
        end
        
        set(gca, 'xtick', 1:10, 'xticklabel', labels);
    end
end

% show distribution of natural images in space of top PCs
figure;
smartscatter(all_pc_proj(:, 1), all_pc_proj(:, 2));
xlabel('PC1');
ylabel('PC2');

%% Get a covariance matrix in the full 30d NI texture space

large_cov = cov(ni_all_ev);

% psychophysical sensitivities \propto 1/distance(texture, random)
% distance(texture, random) = \sum_i C_{ij}
%       (texture_i - random_i) (texture_j - random_j)

%% Converting from binary patches to 30d Chubbified continuous stats

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
        loc_ev = zeros(n_samples_binary, 10*size(chubb_nonlin.pvals, 2));
        for p = 1:n_samples_binary
            crt_patch = crt_patches{j}(:, :, p);
            for k = 1:size(chubb_nonlin.pvals, 2)
                chubbed_patch = applyNonlinearity(crt_patch, chubb_nonlin.pvals(:, k));
                [~, loc_ev(p, 10*(k-1)+1:10*k)] = processBlock(chubbed_patch, inf);
            end
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

%% Converting from ternary patches to 30d Chubbified continuous stats

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
        loc_ev = zeros(n_samples_ternary, 10*size(chubb_nonlin.pvals, 2));
        for p = 1:n_samples_ternary
            crt_patch = crt_patches{j}(:, :, p);
            for k = 1:size(chubb_nonlin.pvals, 2)
                chubbed_patch = applyNonlinearity(crt_patch, chubb_nonlin.pvals(:, k));
                [~, loc_ev(p, 10*(k-1)+1:10*k)] = processBlock(chubbed_patch, inf);
            end
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

%% Diagnose the behavior of Chubbified stats near iid random patches

patch_rnd = ternary_map_patches{1}{1}(:, :, 1); % AB_1_1, [0 1 0]
patch_plus = ternary_map_patches{1}{2}(:, :, 1);
patch_minus = ternary_map_patches{7}{2}(:, :, 1); % AB_1_1, [2/3, -1/3, 2/3]

patches = {patch_rnd, patch_plus, patch_minus};
evs = cell(size(patches));
for i = 1:length(patches)
    crt_ev = zeros(1, 10*size(chubb_nonlin.pvals, 2));
    for j = 1:size(chubb_nonlin.pvals, 2)
        crt_patch_nl = applyNonlinearity(patches{i}, chubb_nonlin.pvals(:, j));
        [~, crt_sub_ev] = processBlock(crt_patch_nl, inf);
        crt_ev(10*(j-1)+1:10*j) = crt_sub_ev;
    end
    evs{i} = crt_ev;
end

dot_norm = @(v, w) dot(v, w)/(norm(v)*norm(w));

diff_ev_plus = evs{2} - evs{1};
diff_ev_minus = evs{3} - evs{1};

%% Plot coordinate axes in PC1/2 space

% calculate PCA for natural image patches
all_pc_coeffs = pca(ni_all_ev);

ni_all_ev_proj_nonmc = ni_all_ev*all_pc_coeffs;

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

% !!! Origin of discrete graylevel patches will not match that from
%     continuous ones. The mapped 30d position of a random continuous patch
%     is related to the means of the Chubb nonlinearities. For a binary
%     patch, only the average between the values of the nonlinearities at 0
%     and at 1 will matter; for a ternary patch, the value at 1/2 will be
%     included in the average. And so on, so that we converge towards the
%     continuous answer in the limit n_grayscale_levels -> infinity.

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
[all_pc_coeffs, all_pc_proj, all_pc_vars, ~, all_pc_explained] = pca(ni_all_ev);

% project binary and ternary thresholds
binary_threshold_locs_projected = cellfun(@(v) ...
    v(:)'*all_pc_coeffs, binary_threshold_locs, 'uniform', false);
ternary_threshold_locs_projected = cellfun(@(v) ...
    v(:)'*all_pc_coeffs, ternary_threshold_locs, 'uniform', false);

% show distribution of natural images in space of top PCs
figure;
ni_all_ev_proj_nonmc = ni_all_ev*all_pc_coeffs;
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

trafo_matrix = sqrtm(large_cov);
[binary_predictions, binary_mapped, binary_optim_details, binary_intersect_details] = ...
    optimizePredictedThresholds(binary_thresholds, binary_interpolated, ...
    trafo_matrix, [0 0.5], 'ignore_mask', strcmp(binary_map_groups, 'A_1'));

[ternary_predictions, ternary_mapped, ternary_optim_details, ternary_intersect_details] = ...
    optimizePredictedThresholds(ternary_thresholds, ternary_interpolated, ...
    trafo_matrix, [0 0.5], 'ignore_mask', strcmp(ternary_map_groups, 'A_1'));

%% 

crt_group_idx = 1;
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 18 10];

crt_means = cell2mat(...
    cellfun(@(m) mean(m, 1), binary_map_stats{crt_group_idx}', 'uniform', false));
crt_low = cell2mat(...
    cellfun(@(m) quantile(m, 0.2, 1), binary_map_stats{crt_group_idx}', 'uniform', false));
crt_high = cell2mat(...
    cellfun(@(m) quantile(m, 0.8, 1), binary_map_stats{crt_group_idx}', 'uniform', false));
for j = 1:size(crt_means, 2)
    ax = axes;
    ax.OuterPosition = [mod(j-1, 6)/6, 4/5 - floor((j-1)/6)/5, 1/6, 1/5];
    
%     crt_coeffs = polyfit(binary_map_locs{crt_group_idx}(:), crt_means(:, j), 3);
    crt_coeffs = polyfit(binary_map_locs{crt_group_idx}(:), crt_means(:, j), 2);
    crt_x = linspace(min(binary_map_locs{crt_group_idx}), max(binary_map_locs{crt_group_idx}), 50);
    crt_fit = crt_coeffs(1)*crt_x.^2 + crt_coeffs(2)*crt_x + crt_coeffs(3);
%     crt_fit = crt_coeffs(1)*crt_x.^3 + crt_coeffs(2)*crt_x.^2 + crt_coeffs(3)*crt_x + crt_coeffs(4);
    
    hold on;
    
    errorbar(binary_map_locs{crt_group_idx}, crt_means(:, j), ...
        crt_means(:, j) - crt_low(:, j), crt_high(:, j) - crt_means(:, j), ...
        'marker', 'none', 'linestyle', 'none', 'color', [0.5 0.5 0.5], ...
        'linewidth', 0.5, 'capsize', 0);
    plot(crt_x, crt_fit);
    scatter(binary_map_locs{crt_group_idx}, crt_means(:, j), '.');
    
%     ylim([-1 1]);
%     ylim([-0.2 0.2]);
end

%% 

crt_group_idx = 145;
fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 18 10];

crt_means = cell2mat(...
    cellfun(@(m) mean(m, 1), ternary_map_stats{crt_group_idx}', 'uniform', false));
crt_low = cell2mat(...
    cellfun(@(m) quantile(m, 0.2, 1), ternary_map_stats{crt_group_idx}', 'uniform', false));
crt_high = cell2mat(...
    cellfun(@(m) quantile(m, 0.8, 1), ternary_map_stats{crt_group_idx}', 'uniform', false));
for j = 1:size(crt_means, 2)
    ax = axes;
    ax.OuterPosition = [mod(j-1, 6)/6, 4/5 - floor((j-1)/6)/5, 1/6, 1/5];
    
%     crt_coeffs = polyfit(ternary_map_locs{crt_group_idx}(:), crt_means(:, j), 3);
    crt_coeffs = polyfit(ternary_map_locs{crt_group_idx}(:), crt_means(:, j), 2);
    crt_x = linspace(min(ternary_map_locs{crt_group_idx}), max(ternary_map_locs{crt_group_idx}), 50);
    crt_fit = crt_coeffs(1)*crt_x.^2 + crt_coeffs(2)*crt_x + crt_coeffs(3);
%     crt_fit = crt_coeffs(1)*crt_x.^3 + crt_coeffs(2)*crt_x.^2 + crt_coeffs(3)*crt_x + crt_coeffs(4);
    
    hold on;
    
    errorbar(ternary_map_locs{crt_group_idx}, crt_means(:, j), ...
        crt_means(:, j) - crt_low(:, j), crt_high(:, j) - crt_means(:, j), ...
        'marker', 'none', 'linestyle', 'none', 'color', [0.5 0.5 0.5], ...
        'linewidth', 0.5, 'capsize', 0);
    plot(crt_x, crt_fit);
    scatter(ternary_map_locs{crt_group_idx}, crt_means(:, j), '.');
    
%     ylim([-1 1]);
%     ylim([-0.2 0.2]);
end