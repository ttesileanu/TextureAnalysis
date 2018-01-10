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
use_mc = true;

for i = 1:length(sub_edirs)
    crt_edir = sub_edirs{i};
    crt_edir_data = data_ext0.ds_merged.avg.edirs.(crt_edir);
    crt_group = crt_edir_data.cgroup_names{1};
    crt_uvecs = crt_edir_data.uvecs;
    crt_ndirs = size(crt_uvecs, 1);
    if ~use_mc
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
    
    groups = [groups repmat({crt_group}, 1, crt_ndirs)]; %#ok<AGROW>
    
    % transform 2d coordinates to 3d probabilities
    crt_axes = 1/3 + crt_uvecs*[[2/3 ; -1/3] [-1/3 ; 2/3] [-1/3 ; -1/3]];
    crt_axes_cell = mat2cell(crt_axes, ones(crt_ndirs, 1), 3);
    
    tex_axes = [tex_axes crt_axes_cell(:)']; %#ok<AGROW>
    directions = [directions repmat({'pos'}, 1, length(crt_uvecs))]; %#ok<AGROW>
    
    thresholds = [thresholds crt_thresh(:)']; %#ok<AGROW>
    stdev_ratios = [stdev_ratios crt_stdev_ratios(:)']; %#ok<AGROW>
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

% restrict correlation order
min_corr_order = 2;
max_corr_order = 2;

corr_mask = cellfun(@(s) length(s)-4 >= min_corr_order && length(s)-4 <= max_corr_order, ...
    groups);
groups = groups(corr_mask);
tex_axes = tex_axes(corr_mask);
directions = directions(corr_mask);
thresholds = thresholds(corr_mask);
stdev_ratios = stdev_ratios(corr_mask);

three_g_to_cont_stats = three_g_to_cont_stats(corr_mask);
all_dot_locs = all_dot_locs(corr_mask);

%% Load natural image stats

natural_stats0 = open(fullfile('save', 'natural_nosky_continuous_with_focus.mat'));
res_choice = 4;

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

%% Map G=3 thresholds to continuous space using a linear approximation

threshold_locations = zeros(length(thresholds), 10);
transformed_directions = zeros(length(thresholds), 10);

for i = 1:length(thresholds)
    % get the derivative of the transformed coordinates at the origin
    % XXX we could do this much more efficiently if we generated more
    %     points closer to the origin

    % now all_dot_locs{i} always starts with 0
    % find where the origin is in dot_locs
    idx0 = 1;
    idx1 = 2;
    dot_step = all_dot_locs{i}(idx1) - all_dot_locs{i}(idx0);
    
    crt_derivative = (mean(three_g_to_cont_stats{i}{idx1}, 1) - ...
                      mean(three_g_to_cont_stats{i}{idx0}, 1)) / dot_step;
                  
    threshold_locations(i, :) = crt_derivative*thresholds(i);
    transformed_directions(i, :) = crt_derivative / norm(crt_derivative);
end

%% Show the thresholds in continuous space

scatterTextureStats({threshold_locations});

%% Predict thresholds from Gaussian-approximated natural statistics

% assume that threshold = variance^trafo_pow
% and that this holds in the principal component coordinate system
trafo_pow = -0.25;

predicted_thresholds0 = zeros(size(thresholds));
trafo_matrix = natural_cov^trafo_pow;
for i = 1:length(thresholds)
    predicted_thresholds0(i) = norm(trafo_matrix*transformed_directions(i, :)')^2;
end

% scale thresholds to match experimental values as closely as possible
mask = isfinite(predicted_thresholds0) & isfinite(thresholds);
% exclude A_1 directions because we equalized natural image patches, so
% that direction is not properly represented by our analysis
mask(strcmp(groups, 'A_1')) = false;
scaling_thresholds = dot(predicted_thresholds0(mask), thresholds(mask)) / ...
    norm(predicted_thresholds0(mask))^2;
predicted_thresholds = predicted_thresholds0*scaling_thresholds;

%% Compare predicted to actual thresholds

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(thresholds) & ~isnan(predicted_thresholds));
% exclude A_1 directions because we equalized natural image patches, so
% that direction is not properly represented by our analysis
mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

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

drawfitline(predicted_thresholds(mask), thresholds(mask), 'showci', false, ...
    'corrtype', 'spearman');
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

safe_print(fullfile('figs', 'G3vsContLinear', ['cont_pred_thresholds' save_tag]));

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

safe_print(fullfile('figs', 'G3vsContLinear', ['cont_pred_thresholds_two_fits' save_tag]));

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

safe_print(fullfile('figs', 'G3vsContLinear', ['cont_ellipses' save_tag]));

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

safe_print(fullfile('figs', 'G3vsContLinear', ['cont_vs_g3ni' save_tag]));

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

safe_print(fullfile('figs', 'G3vsContLinear', ['cont_vs_g3_residuals' save_tag]));