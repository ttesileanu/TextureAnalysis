% Analyzing three grayscale level-data using continuous texture statistics.

%% Load psychophysics data

f = fopen(fullfile('data', 'three_grayscale_thresholds.csv'));
data0 = textscan(f, '%s', 'delimiter', '');
fclose(f);

% process the data a little

% split at delimiters
groups_sub = strsplit(data0{1}{1}, ',', 'CollapseDelimiters', false);
% copy group name over to empty fields
last_type = '';
for i = 1:length(groups_sub)
    if isempty(groups_sub{i})
        groups_sub{i} = last_type;
    else
        % update format to match that required by the texture generation
        % routines
        digit_mask = find(ismember(groups_sub{i}, '0':'9'));
        if ~isempty(digit_mask)
            letter_part = groups_sub{i}(1:digit_mask(1)-1);
            number_part = cell2mat(arrayfun(...
                @(c) ['_' c], groups_sub{i}(digit_mask(1):end), 'uniform', false));
            groups_sub{i} = [letter_part number_part];
        end
        last_type = groups_sub{i};
    end
end
tex_axes_sub0 = strsplit(data0{1}{2}, ',', 'CollapseDelimiters', false);
% copy axes names over to empty fields
last_axis = '';
for i = 1:length(tex_axes_sub0)
    if isempty(tex_axes_sub0{i})
        tex_axes_sub0{i} = last_axis;
    else
        last_axis = tex_axes_sub0{i};
    end
end
tex_axes_sub = cellfun(@(s) arrayfun(@str2double, s(2:end-1)), tex_axes_sub0, 'uniform', false);
directions_sub = strsplit(data0{1}{3}, ',', 'CollapseDelimiters', false);
thresholds_sub = strsplit(data0{1}{4}, ',', 'CollapseDelimiters', false);
stdev_ratios_sub = strsplit(data0{1}{5}, ',', 'CollapseDelimiters', false);
if length(groups_sub) ~= length(tex_axes_sub) || length(groups_sub) ~= length(thresholds_sub) || ...
        length(groups_sub) ~= length(directions_sub) || length(groups_sub) ~= length(stdev_ratios_sub)
    error('Number of elements is not consistent between rows!');
end

% get rid of row names
groups_sub = groups_sub(2:end);
tex_axes_sub = tex_axes_sub(2:end);
directions_sub = directions_sub(2:end);
thresholds_sub = cellfun(@str2double, thresholds_sub(2:end));
stdev_ratios_sub = cellfun(@str2double, stdev_ratios_sub(2:end));

%% Load psychophysics data (extended)

data_ext0 = open(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% select single planes, full data
all_edirs = fieldnames(data_ext0.ds_merged.avg.edirs);
edir_mask = true(size(all_edirs));
for i = 1:length(all_edirs)
    crt_edir = all_edirs{i};
    crt_edir_data = data_ext0.ds_merged.avg.edirs.(crt_edir);
    
    crt_cgroup_names = crt_edir_data.cgroup_names;
    if ~strcmp(crt_cgroup_names{1}, crt_cgroup_names{2})
        % this is not single plane; skip
        edir_mask(i) = false;
        continue;
    end
    
    if ismember([crt_edir 'full'], all_edirs)
        % keep only full
        edir_mask(i) = false;
        continue;
    end

%     % keep only fields not ending in 'full'
%     if length(crt_edir) > 4 && strcmp(crt_edir(end-3:end), 'full')
%         edir_mask(i) = false;
%         continue;
%     end
    
    if ismember(crt_edir, {'A1G', 'A1Ginter'})
        % keep only merge
        edir_mask(i) = false;
        continue;
    end
end
sub_edirs = all_edirs(edir_mask);

% convert data to something similar to the CSV file
groups = {};
tex_axes = {};
directions = {};
thresholds = [];
stdev_ratios = [];
% XXX jwb doesn't have any single-plane measurements
subjects = setdiff(fieldnames(data_ext0.ds_merged), {'avg', 'jwb'});
thresholds_by_subject = cell(1, length(subjects));
stdev_ratios_by_subject = cell(1, length(subjects));
for i = 1:length(sub_edirs)
    crt_edir = sub_edirs{i};
    crt_edir_data = data_ext0.ds_merged.avg.edirs.(crt_edir);
    crt_group = crt_edir_data.cgroup_names{1};
    crt_uvecs = crt_edir_data.uvecs;
    crt_ndirs = size(crt_uvecs, 1);
    
    % get averaged data
    crt_thresh = crt_edir_data.thresh_mags_all;
    crt_thresh_hi = crt_edir_data.thresh_mags_ebhi_all;
    crt_thresh_lo = crt_edir_data.thresh_mags_eblo_all;
    crt_stdev_ratios = crt_thresh_hi./crt_thresh_lo;
    
    thresholds = [thresholds crt_thresh(:)']; %#ok<AGROW>
    stdev_ratios = [stdev_ratios crt_stdev_ratios(:)']; %#ok<AGROW>

    % get per-subject data
    for j = 1:length(subjects)
        thresh_field = ['thresh_mags_' subjects{j}];
        thresh_hi_field = ['thresh_mags_ebhi_' subjects{j}];
        thresh_lo_field = ['thresh_mags_eblo_' subjects{j}];
        
        if all(isfield(crt_edir_data, {thresh_field, thresh_hi_field, thresh_lo_field}))        
            crt_thresh_subj = crt_edir_data.(thresh_field);
            crt_thresh_hi_subj = crt_edir_data.(thresh_hi_field);
            crt_thresh_lo_subj = crt_edir_data.(thresh_lo_field);
            crt_stdev_ratios_subj = crt_thresh_hi_subj./crt_thresh_lo_subj;
        else
            crt_thresh_subj = nan(size(crt_thresh));
            crt_thresh_hi_subj = nan(size(crt_thresh_hi));
            crt_thresh_lo_subj = nan(size(crt_thresh_lo));
            crt_stdev_ratios_subj = nan(size(crt_stdev_ratios));
        end
        
        thresholds_by_subject{j} = [thresholds_by_subject{j} crt_thresh_subj(:)'];
        stdev_ratios_by_subject{j} = [stdev_ratios_by_subject{j} crt_stdev_ratios(:)'];
    end
    
    groups = [groups repmat({crt_group}, 1, crt_ndirs)]; %#ok<AGROW>
    
    % transform 2d coordinates to 3d probabilities
    crt_axes = 1/3 + crt_uvecs*[[2/3 ; -1/3] [-1/3 ; 2/3] [-1/3 ; -1/3]];
    crt_axes_cell = mat2cell(crt_axes, ones(crt_ndirs, 1), 3);
    
    tex_axes = [tex_axes crt_axes_cell(:)']; %#ok<AGROW>
    directions = [directions repmat({'pos'}, 1, length(crt_uvecs))]; %#ok<AGROW>    
end

%% Compare CSV data with the extended data

thresholds_sub_from_ext = zeros(size(thresholds_sub));
stdev_ratio_sub_from_ext = zeros(size(stdev_ratios_sub));
for i = 1:length(thresholds_sub_from_ext)
    crt_group = groups_sub{i};
    crt_axis = tex_axes_sub{i};
    crt_dir = directions_sub{i};
    
    crt_ext_idxs_group = strcmp(groups, crt_group);
    
    % this is by how much we need to multiply extended threshold to obtain
    % value compatible to CSV threshold
    crt_scaling = 1;
    if max(abs(crt_axis - [1 0 0])) < eps
        if strcmp(crt_dir, 'pos')
            crt_axis_ext = [1 0 0];
        else
            crt_axis_ext = [-1/3 2/3 2/3];
        end
    elseif max(abs(crt_axis - [0 1 0])) < eps
        if strcmp(crt_dir, 'pos')
            crt_axis_ext = [0 1 0];
        else
            crt_axis_ext = [2/3 -1/3 2/3];
        end
    elseif max(abs(crt_axis - [0 0 1])) < eps
        crt_scaling = 1/sqrt(2);
        if strcmp(crt_dir, 'pos')
            crt_axis_ext = 1/3 + [-1/3 -1/3 2/3]/sqrt(2);
        else
            crt_axis_ext = 1/3 - [-1/3 -1/3 2/3]/sqrt(2);
        end
    end
    
    crt_ext_idxs_axis = cellfun(@(v) max(abs(v - crt_axis_ext)) < eps, tex_axes);
    
    crt_ext_idxs = find(crt_ext_idxs_group & crt_ext_idxs_axis);
    if isempty(crt_ext_idxs)
        warning(['No match found at index ' int2str(i) '.']);
        break;
    end
    if length(crt_ext_idxs) > 1
        warning(['Ambiguous data found at index ' int2str(i) '.']);
        break;
    end
    
    thresholds_sub_from_ext(i) = thresholds(crt_ext_idxs)*crt_scaling;
    stdev_ratio_sub_from_ext(i) = stdev_ratios(crt_ext_idxs)*crt_scaling;
end

%% Generate patches in the directions used in the psychophysics

% number of patches to generate for each set of coordinates
n_per_dot = 16;

% dot_locs = linspace(-1/2, 1, 16);

% number of dots to generate in each direction
n_dots = 10;

patch_size = 64;

% initialize the texture generator
% overall map setup
btc_dict = btc_define;
%aug_opts = [];
%aug_opts.ifstd = 1;
btc_makemaps_opts = [];
btc_makemaps_opts.nmaps = 1;
btc_makemaps_opts.show = 0;
btc_makemaps_opts.area = [patch_size patch_size];
auxopts = btc_auxopts; % won't need metropolis arguments here

% 3 grayscale levels
mtcs = mtc_define(3);

t0 = tic;
three_g_to_cont_stats = cell(1, length(tex_axes));
all_dot_locs = cell(1, length(tex_axes));
for i = 1:length(tex_axes)
    percentage = (i-1)/length(tex_axes)*100;
    disp(['Working on ' groups{i} ' [' num2str(tex_axes{i}, '%6.2f') ']... (' int2str(round(percentage)) '% done, ' ...
        num2str(toc(t0), '%.2f') ' seconds elapsed)']);
    
    crt_group = groups{i};
    crt_axis = tex_axes{i};
   
    % look at locations dot_locs in this direction
    max_dot_loc = 1/(1 - 3*min(crt_axis));
    crt_dot_locs = linspace(0, max_dot_loc, n_dots);
    all_dot_locs{i} = crt_dot_locs;
    three_g_to_cont_stats{i} = cell(1, n_dots);
    for j = 1:length(crt_dot_locs)
        crt_t = crt_dot_locs(j);
        crt_coords = [1/3 1/3 1/3]*(1 - crt_t) + crt_axis*crt_t;
        
        % set up the texture generator for this particular direction
        cgs_struct = [];
        cgs_struct.(crt_group) = crt_coords;
        augcoords = mtc_augcoords(cgs_struct, mtcs, btc_dict);
        
        % generate a number of patches at each location, and calculate the
        % corresponding continuous stats
        crt_cont_stats = zeros(n_per_dot, 10);
        for k = 1:n_per_dot
            % patches are generated with values 0, 1, 2, so we divide by 2
            % to normalize
            crt_patch = btc_makemaps(augcoords.method{1}, ...
                btc_makemaps_opts, btc_dict, auxopts, []) / 2.0;
            % get continuous stats for the patch
            % XXX shouldn't we equalize first?
            [~, crt_cont_stats(k, :)] = processBlock(crt_patch, inf);
        end
        three_g_to_cont_stats{i}{j} = crt_cont_stats;
    end
end
disp('Done.');

%% Save results

save(fullfile('save', 'three_psycho_via_continuous_extended.mat'), ...
    'directions', 'groups', 'n_per_dot', 'n_dots', 'patch_size', ...
    'tex_axes', 'three_g_to_cont_stats', 'thresholds', 'stdev_ratios', ...
    'all_dot_locs');

%% Load results

load(fullfile('save', 'three_psycho_via_continuous_extended.mat'));

%% Load natural image stats

% natural_stats0 = open(fullfile('save', 'natural_nosky_continuous_with_focus.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_continuous_with_focus_contrastadapt.mat'));
res_choice = 4;

% natural_stats0 = open(fullfile('save', 'natural_nosky_ternarized_continuous_with_focus.mat'));
% res_choice = 4;

%% Preprocessing

restrict_to_focus = true;
% scale_to_binary = true;
scale_to_binary = false;

natural_stats = natural_stats0.res{res_choice};

% scale the continuous stats dimensions to match the variability from a
% binary natural image analysis
if scale_to_binary
    natural_stats_bin0 = open(fullfile('save', 'natural_nosky_binary_with_focus.mat'));
    natural_stats_bin = natural_stats_bin0.res{res_choice};
    ni_scalings = std(natural_stats_bin.ev, [], 1) ./ std(natural_stats.ev, [], 1);
    
    natural_stats.ev = bsxfun(@times, natural_stats.ev, ni_scalings);
    natural_stats.covM = cov(natural_stats.ev);
end

% restrict to in-focus patches, if possible
if restrict_to_focus && isfield(natural_stats, 'focus')
    mask = (natural_stats.focus.clusterIds == natural_stats.focus.focusCluster);
    fields = {'objIds', 'ev', 'pxPerPatch', 'patchLocations', 'patchLocationsOrig', 'imgIds'};
    for i = 1:length(fields)
        natural_stats.(fields{i}) = natural_stats.(fields{i})(mask, :);
    end
    natural_stats.covM = cov(natural_stats.ev);
end

save_tag = '_';
if restrict_to_focus
    save_tag = [save_tag 'focus_'];
else
    save_tag = [save_tag 'all_patches_']; %#ok<*UNRCH>
end
if scale_to_binary
    save_tag = [save_tag 'scaled_'];
else
    save_tag = [save_tag 'not_scaled_'];
end

% save_tag = sprintf('_all_patches_not_scaled_N%d_R%d', natural_stats.options.blockAF, ...
%     natural_stats.options.patchSize(1));
% save_tag = sprintf('_focus_not_scaled_N%d_R%d', natural_stats.options.blockAF, ...
%     natural_stats.options.patchSize(1));
% save_tag = sprintf('_all_patches_scaled_N%d_R%d', natural_stats.options.blockAF, ...
%     natural_stats.options.patchSize(1));
save_tag = [save_tag sprintf('N%d_R%d', natural_stats.options.blockAF, ...
    natural_stats.options.patchSize(1))];

% fit a Gaussian to the data
natural_mu = mean(natural_stats.ev, 1);
natural_cov = cov(natural_stats.ev);
natural_sqrtcov = sqrtm(natural_cov);

%% Map G=3 thresholds to continuous space

threshold_locations = cell(1, length(thresholds));

% use linear interpolation or extrapolation
for i = 1:length(thresholds)
    crt_thresh = thresholds(i);    
    crt_dot_locs = all_dot_locs{i};
    
    % can we interpolate, or do we have to extrapolate?
    k = find(crt_dot_locs(1:end-1) <= crt_thresh & crt_dot_locs(2:end) > crt_thresh, 1);
    if ~isempty(k)
        % yes, we can interpolate
        loc1 = mean(three_g_to_cont_stats{i}{k}, 1);
        loc2 = mean(three_g_to_cont_stats{i}{k+1}, 1);
        
        alpha = (crt_thresh - crt_dot_locs(k)) / (crt_dot_locs(k+1) - crt_dot_locs(k));
        thresh_loc = (1 - alpha)*loc1 + alpha*loc2;
        
        threshold_locations{i} = thresh_loc;
    else
        if crt_thresh < crt_dot_locs(1)
            error('This shouldn''t happen. -- crt_thresh should be >= 0, and crt_dot_locs(1) should be 0.');
        end
        % we need to extrapolate
        loc1 = mean(three_g_to_cont_stats{i}{end-1}, 1);
        loc2 = mean(three_g_to_cont_stats{i}{end}, 1);
        der = (loc2 - loc1) / (crt_dot_locs(end) - crt_dot_locs(end - 1));
        thresh_loc = loc2 + der*(crt_thresh - crt_dot_locs(end));
        threshold_locations{i} = thresh_loc;
    end
end

%% Plot the thresholds in some 10d coordinate planes

% plot the natural image distribution in the 10d continuous texture space,
% with the mapped G=3 thresholds shown on top
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

pairs = {[2 10], [2 3], [3 5], [6 8], [6 10], [4 7]};

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 12 6];

threshold_locations_mat = cell2mat(threshold_locations');

for i = 1:length(pairs)
    ax = axes;
    px = mod(i-1, 3)/3;
    py = floor((i-1) / 3)/2;
    ax.OuterPosition = [px 1/2 - py 1/3 1/2];
    
    ev1 = natural_stats.ev(:, pairs{i}(1));
    ev2 = natural_stats.ev(:, pairs{i}(2));
    smartscatter(ev1, ev2);
    
    hold on;
    plot(threshold_locations_mat(:, pairs{i}(1)), ...
         threshold_locations_mat(:, pairs{i}(2)), '.r');
    
    xlabel(labels{pairs{i}(1)});
    ylabel(labels{pairs{i}(2)});
    
    beautifygraph;
end

preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['thresholds_in_10d' save_tag]), 'png');

%% Optimize match between predicted and measured thresholds

ignore_mask = true(size(thresholds));
% ignore contrast plane because we equalize that out of natural images
ignore_mask(strcmp(groups, 'A_1')) = false;
% % exclude directions for which the psycophysics has infinite error bars
% ignore_mask(isnan(stdev_ratios)) = false;

% range of noise radius^2 in which to search for an optimum
noise2_range = [0 0.03];

optim_options = optimset('display', 'iter', 'tolx', 1e-10);
% calculate the L2 norm of the difference between two vectors while
% ignoring entries that are infinite or NaN in either of them
normdiff_mask = @(a, b, mask) norm((a(mask) - b(mask)) / (sum(mask) - 1));
normdiff_skipinf = @(a, b) normdiff_mask(a, b, ignore_mask & isfinite(a) & isfinite(b));

% % calculate the Spearman correlation between two vectors while ignoring
% % entries that are infinite or NaN in either of them
% spcorr_mask = @(a, b, mask) corr(a(mask), b(mask), 'type', 'spearman');
% spcorr_skipinf = @(a, b) spcorr_mask(a(:), b(:), ignore_mask & isfinite(a) & isfinite(b));

[optimal_size, optimal_norm, ~, optim_details] = fminbnd(@(sz) normdiff_skipinf(...
    predictThreeGThresholds(sz, three_g_to_cont_stats, natural_sqrtcov, all_dot_locs), ...
    thresholds), noise2_range(1), noise2_range(2), optim_options);

% recalculate the predicted thresholds at the optimal noise radius
[predicted_thresholds, predicted_thresh_locs] = predictThreeGThresholds(...
    optimal_size, three_g_to_cont_stats, natural_sqrtcov, all_dot_locs);

%% Make histograms of natural image stats along rays

sigma = 0.1;
natural_stats_hist = cell(size(three_g_to_cont_stats));

n_plots = length(three_g_to_cont_stats);
nx = ceil(sqrt(n_plots));
ny = ceil(n_plots / nx);

aspect = 1.5;
max_x = 19.5;
max_y = 12;

if max_x / nx < aspect * max_y / ny
    fig_x = max_x;
    fig_y = ny * max_x / nx / aspect;
else
    fig_y = max_y;
    fig_x = nx * max_y / ny * aspect;
end

fig = figure;
fig.Units = 'inches';
fig.Position = [0.5 0.5 fig_x fig_y];

for i = 1:length(tex_axes)
    crt_dot_locs = all_dot_locs{i};
    natural_stats_hist{i} = zeros(1, length(crt_dot_locs));
    for j = 1:length(crt_dot_locs)
        mid_point = mean(three_g_to_cont_stats{i}{j}, 1);
        distances = bsxfun(@minus, natural_stats.ev, mid_point);
        counts = exp(-0.5*sum((distances/sigma) .^ 2, 2));
        natural_stats_hist{i}(j) = sum(counts);
    end

    ax = axes;
    ax.Units = 'normalized';
    ax.OuterPosition = [mod(i - 1, nx)/nx ...
        floor((i - 1) / nx)/ny 1/nx 1/ny];
    
    hold on;
    fill([crt_dot_locs(:) ; flipud(crt_dot_locs(:))], ...
        [natural_stats_hist{i}(:) ; zeros(length(crt_dot_locs), 1)], [0.7 0.7 0.7]);
    
    lims = ylim;
    crt_thresh = thresholds(i);
    crt_std_ratio = stdev_ratios(i);
    
    crt_t1 = crt_thresh / crt_std_ratio;
    crt_t2 = crt_thresh * crt_std_ratio;
    
    h_patch = fill([crt_t1 crt_t1 crt_t2 crt_t2], [lims(1) lims(2) lims(2) lims(1)], ...
        [0.7 0.2 0.2]);
    h_patch.FaceAlpha = 0.12;
    
    plot([crt_thresh crt_thresh], [lims(1) lims(2)], 'r--');
    
    xlim([0, 2]);
end

%% Compare predicted to actual thresholds

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(thresholds) & ~isnan(predicted_thresholds));
% exclude A_1 directions because we equalized/contrast-adapted natural
% image patches, so that direction is not properly represented by our
% analysis
mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];

masked_groups = groups(mask);
unique_groups = fliplr(unique(masked_groups));
n_unique_groups = length(unique_groups);
n_prism = size(unique(prism, 'rows'), 1);
prism_colors = prism(n_prism);
n_lines = size(unique(lines, 'rows'), 1);
line_colors = lines(n_lines);
if n_unique_groups <= n_prism
    group_colors = prism_colors(1:n_unique_groups, :);
elseif n_unique_groups <= n_prism + n_lines
    group_colors = [prism_colors ; ...
                    line_colors(1:(n_unique_groups - n_prism), :)];
else
    group_colors = parula(n_unique_groups);
end
threshold_colors = zeros(sum(mask), 3);
for i = 1:n_unique_groups
    for k = 1:3
        threshold_colors(strcmp(masked_groups, unique_groups{i}), k) = group_colors(i, k);
    end
end

hold on;
masked_predicted_thresholds = predicted_thresholds(mask);
masked_thresholds = thresholds(mask);
for i = 1:n_unique_groups
    sub_mask_group = strcmp(masked_groups, unique_groups{i});
    sub_preds = masked_predicted_thresholds(sub_mask_group);
    sub_thresh = masked_thresholds(sub_mask_group);
    
    crt_color = 0.7 + 0.3*group_colors(i, :);
    
    for j = 1:length(sub_preds)
        for k = j+1:length(sub_preds)
            line(sub_preds([j k]), sub_thresh([j k]), 'linewidth', 0.5, ...
                'color', crt_color);
        end
    end
end

thresholds_std_pos = thresholds.*(stdev_ratios - 1);
thresholds_std_neg = thresholds - thresholds./stdev_ratios;
h = errorbar(predicted_thresholds(mask), thresholds(mask), ...
    thresholds_std_neg(mask), thresholds_std_pos(mask), ...
    'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
h.CapSize = 0;
th_h = smartscatter(predicted_thresholds(mask), thresholds(mask), 'color', threshold_colors, ...
    'density', false);
% smartscatter(predicted_thresholds(mask), thresholds(mask), 'color', [1 0 0], ...
%     'density', false);

t_min = min(thresholds(mask));
t_max = max(thresholds(mask));

ylim([0.8*t_min 1.2*t_max]);
axis equal;

% drawfitline(predicted_thresholds(mask), thresholds(mask), 'showci', false, ...
%     'corrtype', 'spearman');
drawfitline(predicted_thresholds(mask), thresholds(mask), 'line', [1 0], ...
    'style', {'--k'}, 'legendloc', 'northwest');
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

for i = 1:length(thresholds)
    if ~mask(i)
        continue;
    end
    text(predicted_thresholds(i)+0.01, thresholds(i), ...
        groups{i}, ...
        'fontsize', 6);        
%         [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']'], ...
end

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsContContrastAdapt', ['cont_pred_thresholds' save_tag]));

%% Compare predicted to actual thresholds, centroids and ellipses, by subject

mask_by_subject = cell(size(subjects));
all_groups = [];
for i = 1:length(subjects)
    % exclude directions for which the psychophysics yielded no thresholds
    mask_by_subject{i} = (~isnan(thresholds_by_subject{i}) & ...
        ~isnan(predicted_thresholds));
    % exclude A_1 directions because we equalized/contrast-adapted natural
    % image patches, so that direction is not properly represented by our
    % analysis
    mask_by_subject{i}(strcmp(groups, 'A_1')) = false;
    % exclude directions for which the psycophysics has infinite error bars
    mask_by_subject{i}(isnan(stdev_ratios_by_subject{i})) = false;
    
    all_groups = [all_groups groups(mask_by_subject{i})]; %#ok<AGROW>
end
unique_groups = fliplr(unique(all_groups));
n_unique_groups = length(unique_groups);
n_prism = size(unique(prism, 'rows'), 1);
prism_colors = prism(n_prism);
n_lines = size(unique(lines, 'rows'), 1);
line_colors = lines(n_lines);
if n_unique_groups <= n_prism
    group_colors = prism_colors(1:n_unique_groups, :);
elseif n_unique_groups <= n_prism + n_lines
    group_colors = [prism_colors ; ...
                    line_colors(1:(n_unique_groups - n_prism), :)];
else
    group_colors = parula(n_unique_groups);
end
group_markers = {{'+', 80}, {'o', 80}, {'*', 80}, {'.', 300}, ...
    {'x', 80}, {'s', 80}, {'d', 80}, {'^', 80}, {'v', 80}, ...
    {'>', 80}, {'<', 80}, {'p', 80}, {'h', 80}};

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 7.5 7];

hold on;

t_min_by_subject = zeros(length(subjects), 1);
t_max_by_subject = zeros(length(subjects), 1);

mean_corrs = zeros(length(subjects), 1);
group_handles = zeros(n_unique_groups, 1);

for k = 1:length(subjects)
    mask = mask_by_subject{k};
    masked_groups = groups(mask);
    threshold_colors = zeros(sum(mask), 3);
    for i = 1:n_unique_groups
        for j = 1:3
            threshold_colors(strcmp(masked_groups, unique_groups{i}), j) = group_colors(i, j);
        end
    end

    masked_predicted_thresholds = predicted_thresholds(mask);
    masked_thresholds = thresholds_by_subject{k}(mask);
    all_mean_preds = zeros(n_unique_groups, 1);
    all_mean_thresh = zeros(n_unique_groups, 1);
    for i = 1:n_unique_groups
        sub_mask_group = strcmp(masked_groups, unique_groups{i});
        sub_preds = masked_predicted_thresholds(sub_mask_group);
        sub_thresh = masked_thresholds(sub_mask_group);
        
        if isempty(sub_preds)
            continue;
        end
        
        crt_color0 = group_colors(i, :);
        crt_color = 0.7 + 0.3*crt_color0;
        
        all_mean_preds(i) = mean(sub_preds);
        all_mean_thresh(i) = mean(sub_thresh);
        
        if sum(sub_mask_group) > 1
            crt_cov = cov([sub_preds(:) sub_thresh(:)]);
            ellipse(all_mean_preds(i), all_mean_thresh(i), crt_cov, 'color', crt_color);
        end
        
        crt_marker = group_markers{mod(i-1, length(group_markers))+1};
        
        crt_h = scatter(all_mean_preds(i), all_mean_thresh(i), crt_marker{2}, crt_marker{1}, ...
            'markerfacecolor', 0.8*crt_color0, ...
            'markeredgecolor', 0.8*crt_color0, ...
            'markerfacealpha', 0.6, ...
            'markeredgealpha', 0.6, ...
            'linewidth', 2);
        if k == 1
            group_handles(i) = crt_h;
        end
    end
    
    corr_mask = (isfinite(all_mean_preds) & isfinite(all_mean_thresh));
    mean_corrs(k) = corr(all_mean_preds(corr_mask), all_mean_thresh(corr_mask));
    
%     for i = 1:length(thresholds_by_subject{k})
%         if ~mask(i)
%             continue;
%         end
%         text(predicted_thresholds(i)+0.01, thresholds_by_subject{k}(i), ...
%             groups{i}, ...
%             'fontsize', 6);
%         %         [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']'], ...
%     end

%     t_min_by_subject(k) = min(masked_thresholds);
%     t_max_by_subject(k) = max(masked_thresholds);

    t_min_by_subject(k) = min(all_mean_thresh(corr_mask));
    t_max_by_subject(k) = max(all_mean_thresh(corr_mask));
end

[~, leg_h] = legend(group_handles, unique_groups, 'location', 'southeast');
for i = 1:n_unique_groups
    crt_h = leg_h(length(leg_h)/2 + i);
    crt_h.Children.MarkerSize = group_markers{i}{2}/10;
    crt_h.Children.LineWidth = 1;
end

% thresholds_std_pos = thresholds.*(stdev_ratios - 1);
% thresholds_std_neg = thresholds - thresholds./stdev_ratios;
% h = errorbar(predicted_thresholds(mask), thresholds(mask), ...
%     thresholds_std_neg(mask), thresholds_std_pos(mask), ...
%     'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
% h.CapSize = 0;
% th_h = smartscatter(predicted_thresholds(mask), thresholds(mask), 'color', threshold_colors, ...
%     'density', false);
% % smartscatter(predicted_thresholds(mask), thresholds(mask), 'color', [1 0 0], ...
% %     'density', false);

t_min = min(t_min_by_subject);
t_max = max(t_max_by_subject);
t_range = [0.9*t_min 1.1*t_max];

plot(t_range, t_range, '--k');
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

axis equal;
xlim(t_range);
ylim(t_range);

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsContContrastAdapt', ['cont_pred_thresholds_means_' save_tag]));

%% Compare predicted to actual thresholds, separate fits by order

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(thresholds) & ~isnan(predicted_thresholds));
% exclude A_1 directions because we equalized natural image patches, so
% that direction is not properly represented by our analysis
mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

low_order_mask = cellfun(@(s) length(s) == 6, groups) & mask;
high_order_mask = cellfun(@(s) length(s) > 6, groups) & mask;

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 10 7];

thresholds_std_pos = thresholds.*(stdev_ratios - 1);
thresholds_std_neg = thresholds - thresholds./stdev_ratios;
h = errorbar(predicted_thresholds(mask), thresholds(mask), ...
    thresholds_std_neg(mask), thresholds_std_pos(mask), ...
    'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
h.CapSize = 0;
hold on;
smartscatter(predicted_thresholds(mask), thresholds(mask), 'color', [1 0 0], ...
    'density', false);

t_min = min(thresholds(mask));
t_max = max(thresholds(mask));

ylim([0.8*t_min 1.2*t_max]);

drawfitline(predicted_thresholds(low_order_mask), thresholds(low_order_mask), ...
    'showci', false, 'corrtype', 'spearman', 'legendloc', 'northwest');
drawfitline(predicted_thresholds(high_order_mask), thresholds(high_order_mask), ...
    'showci', false, 'corrtype', 'spearman', 'legendloc', 'northeast');
drawfitline(predicted_thresholds(mask), thresholds(mask), 'line', [1 0], ...
    'style', {'--k'});
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

for i = 1:length(thresholds)
    if ~mask(i)
        continue;
    end
    text(predicted_thresholds(i)+0.01, thresholds(i), ...
        groups{i}, ...
        'fontsize', 6);
%         [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']'], ...
end

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['cont_pred_thresholds_two_fits' save_tag]));

%% Make some plots showing the mapping between G=3 and continuous directions

% if true, instead of saving some of the scatterplots to file, all of them
% are shown, while waiting for a keypress after each plot
show_all = true;

if ~show_all
    to_draw_stats = [1, 3, 7, 9, 31, 33, 49];
else
    to_draw_stats = 1:length(three_g_to_cont_stats);
end

for i = 1:length(to_draw_stats)
    crt_i = to_draw_stats(i);
    crt_dot_locs = all_dot_locs{i};
    
    scatterTextureStats({three_g_to_cont_stats{crt_i}, crt_dot_locs}, 'alpha', 0.2);
    
    fig = gcf;
    fig.Name = [int2str(i) ' / ' int2str(length(to_draw_stats))];
    ax = axes;
    ax.Units = 'normalized';
    ax.Position = [0 0.85 1 0.15];
    axis off;
    
    text(0.5, 0.5, ['Sweep from t=' num2str(crt_dot_locs(1), '%.2g') ...
        ' (blue) to t=' num2str(crt_dot_locs(end), '%.2g') ' (red) in direction ' ...
        groups{crt_i} '[' arrayfun(@int2str, tex_axes{crt_i}) ']'], ...
        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
        'fontsize', 12, 'fontweight', 'bold');
    
    crt_group = groups{crt_i}(groups{crt_i} ~= '_');
    crt_ax = arrayfun(@int2str, tex_axes{crt_i});
    filename = ['sweep_G3_in_contspace_' crt_group '_' crt_ax];
    if ~show_all
        safe_print(fullfile('figs', filename));
    else     
        waitforbuttonpress;
        close;
    end
end

%% ... same, but show just the means for each location

% if true, instead of saving some of the scatterplots to file, all of them
% are shown, while waiting for a keypress after each plot
show_all = true;

if ~show_all
    to_draw_stats = [1, 3, 7, 9, 31, 33, 49];
else
    to_draw_stats = 1:length(three_g_to_cont_stats);
end

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

pairs = {[1 10], [2 3], [4 5], [6, 7], [8 9]};

for i = 1:length(to_draw_stats)
    crt_i = to_draw_stats(i);
    crt_stats = three_g_to_cont_stats{crt_i};
    crt_means = cell2mat(cellfun(@(m) mean(m, 1), crt_stats', 'uniform', false));
    crt_stds = cell2mat(cellfun(@(m) std(m, [], 1), crt_stats', 'uniform', false));
    
    crt_dot_locs = all_dot_locs{i};
    
    colors0 = redblue(255);
    colors0 = max(min(colors0 - (sum(colors0, 2) - 1)/2, 1), 0);
    cidxs = min(max(1 + round(0.5*(1 + crt_dot_locs)*(size(colors0, 1) - 1)), 1), size(colors0, 1));
    colors = colors0(cidxs, :);
    
    fig = figure;
    fig.Units = 'inches';
    fig.Name = [int2str(i) ' / ' int2str(length(to_draw_stats))];
    oldPos = fig.Position;
    centroid = oldPos(1:2) + 0.5*oldPos(3:4);
    figSize = [10 1.6];
    fig.Position = [centroid - 0.5*figSize figSize];
    
    ax = cell(1, length(pairs));
    for j = 1:length(pairs)
        idx1 = pairs{j}(1);
        idx2 = pairs{j}(2);
        
        ax{j} = axes;
        ax{j}.Units = 'normalized';
        ax{j}.OuterPosition = [(j-1)/length(pairs) 0 1/length(pairs) 0.9];

        smartscatter(crt_means(:, idx1), crt_means(:, idx2), 'color', colors, 'density', false);
        xlabel(labels{idx1});
        ylabel(labels{idx2});
        beautifygraph;
        
        xlim([-1 1]);
        ylim([-1 1]);
    end
    
    mean_bottom = mean(cellfun(@(a) a.Position(2), ax));
    mean_height = mean(cellfun(@(a) a.Position(4), ax));
    for j = 1:length(ax)
        ax{j}.Position(2) = mean_bottom;
        ax{j}.Position(4) = mean_height;
    end
    
    preparegraph;
    
    ax = axes;
    ax.Units = 'normalized';
    ax.Position = [0 0.85 1 0.15];
    axis off;
    
    text(0.5, 0.5, ['Sweep from t=' num2str(crt_dot_locs(1), '%.2g') ...
        ' (blue) to t=' num2str(crt_dot_locs(end), '%.2g') ' (red) in direction ' ...
        groups{crt_i} '[' arrayfun(@int2str, tex_axes{crt_i}) ']'], ...
        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
        'fontsize', 12, 'fontweight', 'bold');
    
    crt_group = groups{crt_i}(groups{crt_i} ~= '_');
    crt_ax = arrayfun(@int2str, tex_axes{crt_i});
%    filename = ['sweep_G3_in_contspace_' crt_group '_' crt_ax];
    if ~show_all
%        safe_print(fullfile('figs', filename));
    else     
        waitforbuttonpress;
        close;
    end
end

%% Make plots for each texture group

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 18 12];

unique_groups = fliplr(unique(groups));

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

for i = 1:length(unique_groups)
    ax = axes;
    px = mod(i-1, 4);
    py = floor((i-1)/4);
    ax.Units = 'normalized';
    ax.OuterPosition = [px/4, 1 - (py+1)/3, 1/4, 1/3];
    
    hold on;
    
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
    
    % draw measured and predicted thresholds
    mask_idxs = find(strcmp(groups, unique_groups{i}));
    measured_points = [];
    predicted_points = [];
    h_meas = [];
    h_pred = [];
    for j = 1:length(mask_idxs)
        crt_idx = mask_idxs(j);
        crt_axis = tex_axes{crt_idx};
        
        % measured threshold
        crt_threshold = thresholds(crt_idx);
        crt_thresh_pos = [1/3 1/3 1/3]*(1 - crt_threshold) + crt_axis*crt_threshold;
                
        % predicted threshold
        crt_predicted_threshold = predicted_thresholds(crt_idx);
        crt_pred_thresh_pos = [1/3 1/3 1/3]*(1 - crt_predicted_threshold) + ...
            crt_axis*crt_predicted_threshold;
        
        % this is a litle roundabout, but easier conceptually
        crt_t1 = (3*crt_thresh_pos(2) - 1)/2;
        crt_t2 = (crt_thresh_pos(3) - crt_thresh_pos(1)) * max_t2;
        
        crt_pred_t1 = (3*crt_pred_thresh_pos(2) - 1)/2;
        crt_pred_t2 = (crt_pred_thresh_pos(3) - crt_pred_thresh_pos(1)) * max_t2;
        
        % store points for fitting ellipse
        if ~any(isnan([crt_t1 crt_t2]))
            measured_points = [measured_points ; [crt_t1 crt_t2]]; %#ok<AGROW>
        end
        if ~any(isnan([crt_pred_t1 crt_pred_t2]))
            predicted_points = [predicted_points ; [crt_pred_t1 crt_pred_t2]]; %#ok<AGROW>
        end
        
        % plot
        h_meas = [h_meas plot(crt_t1, crt_t2, 'kx', 'linewidth', 1)]; %#ok<AGROW>
        h_pred = [h_pred plot(crt_pred_t1, crt_pred_t2, 'r.')]; %#ok<AGROW>
    end
    
    legend([h_meas(1) h_pred(1)], 'measurements', 'predictions');
    
    % show ellipses
    if size(measured_points, 1) > 2
        measured_M = fit_ellipse(measured_points);
        [crt_V, crt_D] = eig(measured_M);
        crt_D = diag(crt_D);
        crt_D(crt_D < 0 & crt_D > -1e-4) = 0;
        if all(crt_D >= 0)
            crt_y = [1/sqrt(crt_D(1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2))*sin(angle_range(:))]';
            crt_x = crt_V*crt_y;
            plot(crt_x(1, :), crt_x(2, :), 'color', [0.5 0.5 0.5]);
        end
    end
    if size(predicted_points, 1) > 2
        predicted_M = fit_ellipse(predicted_points);
        [crt_V, crt_D] = eig(predicted_M);
        crt_D = diag(crt_D);
        crt_D(crt_D < 0 & crt_D > -1e-4) = 0;
        if all(diag(crt_D) >= 0)
            crt_y = [1/sqrt(crt_D(1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2))*sin(angle_range(:))]';
            crt_x = crt_V*crt_y;
            plot(crt_x(1, :), crt_x(2, :), 'color', [1 0.6 0.6]);
        end
    end
    
    % arrange plot sizes and labels
    axis equal;
    
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5]);
    
    xlim([-2 2]);
    ylim([-2 2]);
    
    title(unique_groups{i});
    
%    beautifygraph;
    set(gca, 'box', 'on', 'xminortick', 'on', 'yminortick', 'on', 'linewidth', 1);
end

preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['cont_ellipses' save_tag]));

%% Check how much thresholds change with noise radius

sizes = linspace(noise2_range(2)/150, noise2_range(2), 100);

thresholds_vs_size = arrayfun(@(sz) ...
    predictThreeGThresholds(sz, three_g_to_cont_stats, natural_sqrtcov, all_dot_locs), ...
    sizes, 'uniform', false);

% some useful definitions
ignore_mask = true(size(thresholds));
ignore_mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% ignore_mask(isnan(stdev_ratios)) = false;

normdiff_mask = @(a, b, mask) norm((a(mask) - b(mask)) / (sum(mask) - 1));
normdiff_skipinf = @(a, b) normdiff_mask(a, b, ignore_mask & isfinite(a) & isfinite(b));
spcorr_mask = @(a, b, mask) corr(a(mask), b(mask), 'type', 'spearman');
spcorr_skipinf = @(a, b) spcorr_mask(a(:), b(:), ignore_mask & isfinite(a) & isfinite(b));
corr_mask = @(a, b, mask) corr(a(mask), b(mask));
corr_skipinf = @(a, b) corr_mask(a(:), b(:), ignore_mask & isfinite(a) & isfinite(b));

% how similar are predicted to measured thresholds, by RMS and correlation
norm_vs_size = cellfun(@(pred) normdiff_skipinf(pred, thresholds), thresholds_vs_size);
spcorr_vs_size = cellfun(@(pred) spcorr_skipinf(pred, thresholds), thresholds_vs_size);
corr_vs_size = cellfun(@(pred) corr_skipinf(pred, thresholds), thresholds_vs_size);

% plot
figure;
hold on;
yyaxis left;
plot(sizes, spcorr_vs_size, 'b');
plot(sizes, corr_vs_size, 'b--');
ylim([0.5 0.9]);
ylabel('Correlation to experiment');

yyaxis right;
plot(sizes, norm_vs_size, 'r');
ylabel('RMS from experiment');

legend({'Rank corr', 'Linear corr', 'RMS'}, 'location' ,'southeast');

xlabel('Noise radius squared r^2');
title('All directions included');

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['noiserad_dependence' save_tag]));

%% Dependence on noise radius focusing on only low-order planes

% keep only second-order planes
low_order_mask = cellfun(@(s) length(s) == 6, groups);

% how similar are predicted to measured thresholds, by RMS and correlation
norm_vs_size_low_order = cellfun(@(pred) normdiff_mask(pred, thresholds, ...
    low_order_mask & isfinite(pred) & isfinite(thresholds)), thresholds_vs_size);
spcorr_vs_size_low_order = cellfun(@(pred) spcorr_mask(pred(:), thresholds(:), ...
    low_order_mask & isfinite(pred) & isfinite(thresholds)), thresholds_vs_size);
corr_vs_size_low_order = cellfun(@(pred) corr_mask(pred(:), thresholds(:), ...
    low_order_mask & isfinite(pred) & isfinite(thresholds)), thresholds_vs_size);

% plot
figure;
hold on;
yyaxis left;
ylim([0.5 0.9]);
plot(sizes, spcorr_vs_size_low_order, 'b');
plot(sizes, corr_vs_size_low_order, 'b--');
ylabel('Correlation to experiment');

yyaxis right;
plot(sizes, norm_vs_size_low_order, 'r');
ylabel('RMS from experiment');

legend({'Rank corr', 'Linear corr', 'RMS'}, 'location' ,'southeast');

xlabel('Noise radius squared r^2');
title('Focusing on second-order correlations');

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['noiserad_dependence_low_order' save_tag]));

%% Get optimal based on low-order

ignore_mask = true(size(thresholds));
ignore_mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% ignore_mask(isnan(stdev_ratios)) = false;

optim_options = optimset('display', 'iter', 'tolx', 1e-10);
normdiff_low_order = @(pred, threholds) normdiff_mask(pred, thresholds, ...
    low_order_mask & isfinite(pred) & isfinite(thresholds));
optimal_size_low_order = fminbnd(@(sz) normdiff_low_order(...
    predictThreeGThresholds(sz, three_g_to_cont_stats, natural_sqrtcov, all_dot_locs), ...
    thresholds), noise2_range(1), noise2_range(2), optim_options);

[predicted_thresholds_low_order, predicted_thresh_locs_low_order] = ...
    predictThreeGThresholds(optimal_size_low_order, three_g_to_cont_stats, ...
    natural_sqrtcov, all_dot_locs);

% how different are they from the other thresholds we obtained?
figure;
scatterfit(predicted_thresholds(ignore_mask), predicted_thresholds_low_order(ignore_mask));
xlabel('Predicted thresholds, optimum based on all data')
ylabel('Predicted thresholds, optimum based on second-order planes')

beautifygraph;
preparegraph;

%% Make plots again

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 18 12];

unique_groups = fliplr(unique(groups));

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

for i = 1:length(unique_groups)
    ax = axes;
    px = mod(i-1, 4);
    py = floor((i-1)/4);
    ax.Units = 'normalized';
    ax.OuterPosition = [px/4, 1 - (py+1)/3, 1/4, 1/3];
%    subplot(3, 4, i);
    
    hold on;
    
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
    text(-1.1, max_t2+0.05, '[0,0,1]', 'fontsize', 12);
    text(-1.1, -max_t2-0.01, '[1,0,0]', 'fontsize', 12);
    
    % draw measured and predicted thresholds
    mask_idxs = find(strcmp(groups, unique_groups{i}));
    measured_points = [];
    predicted_points = [];
    h_meas = [];
    h_pred = [];
    for j = 1:length(mask_idxs)
        crt_idx = mask_idxs(j);
        crt_axis = tex_axes{crt_idx};
        crt_direction = directions{crt_idx};
        
        % measured threshold
        crt_threshold = (-1)^strcmp(crt_direction, 'neg')*thresholds(crt_idx);
        crt_thresh_pos = [1/3 1/3 1/3]*(1 - crt_threshold) + crt_axis*crt_threshold;
                
        % predicted threshold
        crt_predicted_threshold = (-1)^strcmp(crt_direction, 'neg')*predicted_thresholds_low_order(crt_idx);
        crt_pred_thresh_pos = [1/3 1/3 1/3]*(1 - crt_predicted_threshold) + ...
            crt_axis*crt_predicted_threshold;
        
        % this is a litle roundabout, but easier conceptually
        crt_t1 = (3*crt_thresh_pos(2) - 1)/2;
        crt_t2 = (crt_thresh_pos(3) - crt_thresh_pos(1)) * max_t2;
        
        crt_pred_t1 = (3*crt_pred_thresh_pos(2) - 1)/2;
        crt_pred_t2 = (crt_pred_thresh_pos(3) - crt_pred_thresh_pos(1)) * max_t2;
        
        % store points for fitting ellipse
        if ~any(isnan([crt_t1 crt_t2]))
            measured_points = [measured_points ; [crt_t1 crt_t2]]; %#ok<AGROW>
        end
        if ~any(isnan([crt_pred_t1 crt_pred_t2]))
            predicted_points = [predicted_points ; [crt_pred_t1 crt_pred_t2]]; %#ok<AGROW>
        end
        
        % plot
        h_meas = [h_meas plot(crt_t1, crt_t2, 'kx', 'linewidth', 1)]; %#ok<AGROW>
        h_pred = [h_pred plot(crt_pred_t1, crt_pred_t2, 'r.')]; %#ok<AGROW>
    end
    
    legend([h_meas(1) h_pred(1)], 'measurements', 'predictions');
    
    % show ellipses
    if size(measured_points, 1) > 2
        measured_M = fit_ellipse(measured_points);
        [crt_V, crt_D] = eig(measured_M);
        crt_D = diag(crt_D);
        crt_D(crt_D < 0 & crt_D > -1e-4) = 0;
        if all(crt_D >= 0)
            crt_y = [1/sqrt(crt_D(1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2))*sin(angle_range(:))]';
            crt_x = crt_V*crt_y;
            plot(crt_x(1, :), crt_x(2, :), 'color', [0.5 0.5 0.5]);
        end
    end
    if size(predicted_points, 1) > 2
        predicted_M = fit_ellipse(predicted_points);
        [crt_V, crt_D] = eig(predicted_M);
        crt_D = diag(crt_D);
        crt_D(crt_D < 0 & crt_D > -1e-4) = 0;
        if all(diag(crt_D) >= 0)
            crt_y = [1/sqrt(crt_D(1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2))*sin(angle_range(:))]';
            crt_x = crt_V*crt_y;
            plot(crt_x(1, :), crt_x(2, :), 'color', [1 0.6 0.6]);
        end
    end
    
    % arrange plot sizes and labels
    axis equal;
    
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5]);
    
    xlim([-2 2]);
    ylim([-2 2]);
    
    title(unique_groups{i});   
    
%    beautifygraph;
    set(gca, 'box', 'on', 'xminortick', 'on', 'yminortick', 'on', 'linewidth', 1);
end

preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['cont_ellipses_low_order' save_tag]));

%% Comparison to John's data

john = open(fullfile('/Users/ttesileanu/Dropbox/Textures/GrayLevelNIStats/RawNIG3Stats', ...
    'NaturalImageStatsNSL_G3_N2_PC1_32_32.mat'));
john.names = {'A_1', 'AB_1_1', 'AB_1_2', 'AC_1_1', 'AC_1_2', 'BC_1_1', ...
    'BC_1_2', 'AD_1_1', 'AD_1_2', 'ABC_1_1_1', 'ABC_1_2_2', 'ABC_1_2_1', ...
    'ABC_1_1_2', 'ABD_1_1_1', 'ABD_1_2_2', 'ABD_1_2_1', 'ABD_1_1_2', ...
    'ACD_1_1_1', 'ACD_1_2_2', 'ACD_1_2_1', 'ACD_1_1_2', 'BCD_1_1_1', ...
    'BCD_1_2_2', 'BCD_1_2_1', 'BCD_1_1_2', 'ABCD_1_1_1_1', 'ABCD_1_2_2_2', ...
    'ABCD_1_2_1_1', 'ABCD_1_1_2_2', 'ABCD_1_1_2_1', 'ABCD_1_2_1_2', ...
    'ABCD_1_2_2_1', 'ABCD_1_1_1_2'};
john.all_variances = std(john.stats, [], 3);

john_variances = zeros(size(thresholds));
for i = 1:length(john_variances)
    crt_name = groups{i};
    crt_name_idx = find(strcmp(john.names, crt_name));
    if isempty(crt_name_idx)
        error('Group name not found in John''s data.');
    end
    
    crt_ax = tex_axes{i};
    crt_proj = squeeze(sum(bsxfun(@times, john.stats(crt_name_idx, :, :), crt_ax), 2)) / ...
        norm(crt_ax)^2;
    crt_var = std(crt_proj);
    
    if crt_var < 1e-6
        % XXX a hack
        crt_var = 0;
    end
    
    john_variances(i) = crt_var;
end
john_thresholds = 1./sqrt(john_variances);

%% Show how John's thresholds compare to mine

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 10 3];

ignore_mask = true(size(thresholds));
% ignore contrast plane because we equalize that out of natural images
ignore_mask(strcmp(groups, 'A_1')) = false;

ax = axes;
ax.OuterPosition = [0 0 1/3 1];
scatterfit(thresholds(ignore_mask), predicted_thresholds(ignore_mask), ...
    'scatteropts', {'color', [1 0 0]});
xlabel('Measured thresholds');
ylabel('Predicted thresholds');
title('Continuous stats vs. experiment');

beautifygraph;

ax = axes;
ax.OuterPosition = [1/3 0 1/3 1];
scatterfit(thresholds, john_thresholds, 'scatteropts', {'color', [1 0 0]});
xlabel('Measured thresholds');
ylabel('Predicted thresholds');
title('G=3 stats vs. experiment');

beautifygraph;

ax = axes;
ax.OuterPosition = [2/3 0 1/3 1];
scatterfit(john_thresholds, predicted_thresholds, 'scatteropts', {'color', [1 0 0]})
xlabel('G=3 stats predictions');
ylabel('Continuous stats predictions');
title('Continuous stats vs. G=3 stats');

beautifygraph;

preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['cont_vs_g3ni' save_tag]));

%% Show how residuals between theory&experiment compare between John and me

ignore_mask = true(size(thresholds));
% ignore contrast plane because we equalize that out of natural images
ignore_mask(strcmp(groups, 'A_1')) = false;

[~, my_stats] = drawfitline(thresholds(ignore_mask), predicted_thresholds(ignore_mask), ...
    'nodraw', true);
[~, john_stats] = drawfitline(thresholds, john_thresholds, 'nodraw', true);

figure;
scatterfit(john_stats.residuals, my_stats.residuals, 'scatteropts', {'color', [1 0 0]});
xlabel('Residuals continuous predictions');
ylabel('Residuals G=3 predictions');

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'G3vsCont', ['cont_vs_g3_residuals' save_tag]));

%% notes

% at ellipse level:
%  * second order planes match very well between G=3 and
%    continuous stats (the ones I calculated, which are the ones for which
%    there are measurements)
%  * third order planes don't match
%  * fourth order planes match in orientation, but the G=3 thresholds are
%    much smaller than the continuous ones