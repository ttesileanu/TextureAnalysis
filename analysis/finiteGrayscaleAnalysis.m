% analyze the G=3 data using the G=3 stats

%% Load ternary stats

% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog_contrastadapt_with_focus.mat'));
natural_stats = natural_stats0.res{1};

restrict_to_focus = true;

% restrict to in-focus patches, if possible
if restrict_to_focus && isfield(natural_stats, 'focus')
    mask = (natural_stats.focus.clusterIds == natural_stats.focus.focusCluster);
    fields = {'objIds', 'ev', 'pxPerPatch', 'patchLocations', 'patchLocationsOrig', 'imgIds'};
    for i = 1:length(fields)
        natural_stats.(fields{i}) = natural_stats.(fields{i})(mask, :);
    end
    natural_stats.covM = cov(natural_stats.ev);
end

% expand natural stats to contain all 99 probabilities
ni_ev_full = reshape(natural_stats.ev, [], 2, 33);
ni_ev_full(:, 3, :) = 1 - sum(ni_ev_full, 2);
ni_ev_full = reshape(ni_ev_full, [], 99);

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
% XXX jwb doesn't have any single-plane measurements
ternary_subjects = setdiff(fieldnames(ternary_data_ext0.ds_merged), ...
    {'avg', 'jwb'});
ternary_thresholds_by_subject = cell(1, length(ternary_subjects));
ternary_stdev_ratios_by_subject = cell(1, length(ternary_subjects));

for i = 1:length(ternary_sub_edirs)
    crt_edir = ternary_sub_edirs{i};
    crt_edir_data = ternary_data_ext0.ds_merged.avg.edirs.(crt_edir);
    crt_group = crt_edir_data.cgroup_names{1};
    crt_uvecs = crt_edir_data.uvecs;
    crt_ndirs = size(crt_uvecs, 1);
    
    % get averaged data
    crt_thresh = crt_edir_data.thresh_mags_all;
    crt_thresh_hi = crt_edir_data.thresh_mags_ebhi_all;
    crt_thresh_lo = crt_edir_data.thresh_mags_eblo_all;
    crt_stdev_ratios = crt_thresh_hi./crt_thresh_lo;
    
    ternary_thresholds = [ternary_thresholds crt_thresh(:)']; %#ok<AGROW>
    ternary_stdev_ratios = [ternary_stdev_ratios crt_stdev_ratios(:)']; %#ok<AGROW>
    
    % get per-subject data
    for j = 1:length(ternary_subjects)
        thresh_field = ['thresh_mags_' ternary_subjects{j}];
        thresh_hi_field = ['thresh_mags_ebhi_' ternary_subjects{j}];
        thresh_lo_field = ['thresh_mags_eblo_' ternary_subjects{j}];
        
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
        
        ternary_thresholds_by_subject{j} = ...
            [ternary_thresholds_by_subject{j} crt_thresh_subj(:)'];
        ternary_stdev_ratios_by_subject{j} = ...
            [ternary_stdev_ratios_by_subject{j} crt_stdev_ratios_subj(:)'];
    end
    
    ternary_groups = [ternary_groups repmat({crt_group}, 1, crt_ndirs)]; %#ok<AGROW>
    
    % transform 2d coordinates to 3d probabilities
    crt_axes = 1/3 + crt_uvecs*[[2/3 ; -1/3] [-1/3 ; 2/3] [-1/3 ; -1/3]];
    crt_axes_cell = mat2cell(crt_axes, ones(crt_ndirs, 1), 3);
    
    ternary_axes = [ternary_axes crt_axes_cell(:)']; %#ok<AGROW>
end

% XXX
% XXX
% weird correction
% ternary_thresholds = ternary_thresholds ./ cellfun(@(s) 1./norm(s)^2, ternary_axes);
% XXX
% XXX

%% Show natural stats in PC1/2 plane

[pca_coeffs, pca_score, ~, ~, pca_explained] = pca(ni_ev_full);

smartscatter(pca_score(:, 1), pca_score(:, 2));
xlabel('PC1');
ylabel('PC2');

beautifygraph;
preparegraph;

%% Get threshold predictions

ni_cov = cov(ni_ev_full);
% scaling_mat = sqrtm(ni_cov);
scaling_mat = sqrtm(ni_cov + 1e-10*eye(size(ni_cov)));
% scaling_mat = ni_cov;

% get mapping from the 99 probabilities to texture groups
mtc = processBlock('mtc', 3);
unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
groups = repelem(unique_groups, 3);

predictions = zeros(size(ternary_groups));
for i = 1:length(ternary_groups)
    crt_axis = zeros(size(scaling_mat, 1), 1);
    crt_mask = strcmp(groups, ternary_groups{i});
    if sum(crt_mask) ~= 3
        error('Some group appears something different from 3 times.');
    end
    crt_axis(crt_mask) = ternary_axes{i} - 1/3;
%     predictions(i) = sum(crt_axis.^2)/(crt_axis'*scaling_mat*crt_axis);
    predictions(i) = 1/sqrt(crt_axis'*scaling_mat*crt_axis);
end

%% Compare predicted to actual thresholds

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(ternary_thresholds) & ~isnan(predictions));
% exclude A_1 directions because we equalized/contrast-adapted natural
% image patches, so that direction is not properly represented by our
% analysis
mask(strcmp(ternary_groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

% find optimal scaling of predictions
% S = sum((a*pred - thresh).^2)
% need argmin(a) S
% dS/da = sum((a*pred - thresh) .* pred)
%       = a*dot(pred, pred) - dot(thresh, pred) != 0
% a = dot(thresh, pred) / dot(pred, pred)

axis_scaling = cellfun(@(crt_axis) sqrt(3*norm(crt_axis)^2 - 1)/sqrt(2), ternary_axes);

pred_scaling = dot(ternary_thresholds(mask), axis_scaling(mask).*predictions(mask)) / ...
    sum(predictions(mask).^2);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];

masked_groups = ternary_groups(mask);
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
thresholds_std_pos = ternary_thresholds.*(ternary_stdev_ratios - 1);
thresholds_std_neg = ternary_thresholds - ternary_thresholds./ternary_stdev_ratios;
h = errorbar(pred_scaling*axis_scaling(mask).*predictions(mask), ...
    axis_scaling(mask).*ternary_thresholds(mask), ...
    axis_scaling(mask).*thresholds_std_neg(mask), ...
    axis_scaling(mask).*thresholds_std_pos(mask), ...
    'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
h.CapSize = 0;

masked_predictions = predictions(mask);
masked_thresholds = ternary_thresholds(mask);
masked_axis_scaling = axis_scaling(mask);
for i = 1:n_unique_groups
    sub_mask_group = strcmp(masked_groups, unique_groups{i});
    sub_preds = masked_axis_scaling(sub_mask_group).*masked_predictions(sub_mask_group);
    sub_thresh = masked_axis_scaling(sub_mask_group).*masked_thresholds(sub_mask_group);
    
    crt_color = 0.5 + 0.5*group_colors(i, :);
    
    for j = 1:length(sub_preds)
        for k = j+1:length(sub_preds)
            line(pred_scaling*sub_preds([j k]), sub_thresh([j k]), 'linewidth', 0.5, ...
                'color', crt_color);
        end
    end
end

th_h = smartscatter(pred_scaling*axis_scaling(mask).*predictions(mask), ...
    axis_scaling(mask).*ternary_thresholds(mask), 'color', threshold_colors, ...
    'density', false);

t_min = min(axis_scaling(mask).*ternary_thresholds(mask));
t_max = max(axis_scaling(mask).*ternary_thresholds(mask));

ylim([0.8*t_min 1.2*t_max]);
axis equal;

% drawfitline(predicted_thresholds(mask), thresholds(mask), 'showci', false, ...
%     'corrtype', 'spearman');
drawfitline(pred_scaling*axis_scaling(mask).*predictions(mask), ...
    axis_scaling(mask).*ternary_thresholds(mask), 'line', [1 0], ...
    'style', {'--k'}, 'legendloc', 'northwest');
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

for i = 1:length(ternary_thresholds)
    if ~mask(i)
        continue;
    end
    
    text(pred_scaling*axis_scaling(i).*predictions(i)+0.01, ...
        axis_scaling(i).*ternary_thresholds(i), ternary_groups{i}, ...
        'fontsize', 6);        
%         [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']'], ...
end

beautifygraph;
preparegraph;

%% Compare predicted to actual thresholds, centroids and ellipses, by subject

axis_scaling = cellfun(@(crt_axis) sqrt(3*norm(crt_axis)^2 - 1)/sqrt(2), ternary_axes);

pred_scaling = dot(ternary_thresholds(mask), axis_scaling(mask).*predictions(mask)) / ...
    sum(predictions(mask).^2);

mask_by_subject = cell(size(ternary_subjects));
all_groups = [];
for i = 1:length(ternary_subjects)
    % exclude directions for which the psychophysics yielded no thresholds
    mask_by_subject{i} = (~isnan(ternary_thresholds_by_subject{i}) & ...
        ~isnan(predictions));
    % exclude A_1 directions because we equalized/contrast-adapted natural
    % image patches, so that direction is not properly represented by our
    % analysis
    mask_by_subject{i}(strcmp(ternary_groups, 'A_1')) = false;
    % exclude directions for which the psycophysics has infinite error bars
    mask_by_subject{i}(isnan(ternary_stdev_ratios_by_subject{i})) = false;
%     mask_by_subject{i}(cellfun(@length, ternary_groups) == 12) = false;
    
    all_groups = [all_groups ternary_groups(mask_by_subject{i})]; %#ok<AGROW>
end
unique_groups = fliplr(unique(all_groups));
n_unique_groups = length(unique_groups);
n_prism = size(unique(prism, 'rows'), 1);
prism_colors = prism(n_prism);
n_lines = size(unique(lines, 'rows'), 1);
line_colors = lines(n_lines);
% if n_unique_groups <= n_prism
%     group_colors = prism_colors(1:n_unique_groups, :);
% elseif n_unique_groups <= n_prism + n_lines
%     group_colors = [prism_colors ; ...
%                     line_colors(1:(n_unique_groups - n_prism), :)];
% else
%     group_colors = parula(n_unique_groups);
% end
% group_colors = parula(n_unique_groups);
group_colors = jet(n_unique_groups);
group_markers = {{'+', 60}, {'o', 60}, {'*', 60}, {'.', 250}, ...
    {'x', 60}, {'s', 60}, {'d', 60}, {'^', 60}, {'v', 60}, ...
    {'>', 60}, {'<', 60}, {'p', 60}, {'h', 60}};

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5 4.5];

hold on;

t_min_by_subject = zeros(length(ternary_subjects), 1);
t_max_by_subject = zeros(length(ternary_subjects), 1);

mean_corrs = zeros(length(ternary_subjects), 1);
group_handles = zeros(n_unique_groups, 1);

for k = 1:length(ternary_subjects)    
    mask = mask_by_subject{k};
    masked_groups = ternary_groups(mask);
    threshold_colors = zeros(sum(mask), 3);
    for i = 1:n_unique_groups
        for j = 1:3
            threshold_colors(strcmp(masked_groups, unique_groups{i}), j) = group_colors(i, j);
        end
    end

    masked_predicted_thresholds = pred_scaling * axis_scaling(mask) .* ...
        predictions(mask);
    masked_thresholds = ternary_thresholds_by_subject{k}(mask);
    all_mean_preds = nan(n_unique_groups, 1);
    all_mean_thresh = nan(n_unique_groups, 1);
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
%             ellipse(all_mean_preds(i), all_mean_thresh(i), crt_cov, 'color', crt_color);
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
xlabel('Predicted average thresholds');
ylabel('Measured average thresholds');

axis equal;
xlim(t_range);
ylim(t_range);

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'finite_g3_pred_means'), 'png');

%% Compare predicted to actual thresholds, just averages per plane

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(ternary_thresholds) & ~isnan(predictions));
% exclude A_1 directions because we equalized/contrast-adapted natural
% image patches, so that direction is not properly represented by our
% analysis
mask(strcmp(ternary_groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

% find optimal scaling of predictions
% S = sum((a*pred - thresh).^2)
% need argmin(a) S
% dS/da = sum((a*pred - thresh) .* pred)
%       = a*dot(pred, pred) - dot(thresh, pred) != 0
% a = dot(thresh, pred) / dot(pred, pred)

axis_scaling = cellfun(@(crt_axis) sqrt(3*norm(crt_axis)^2 - 1)/sqrt(2), ternary_axes);

pred_scaling = dot(ternary_thresholds(mask), axis_scaling(mask).*predictions(mask)) / ...
    sum(predictions(mask).^2);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];

masked_groups = ternary_groups(mask);
unique_groups = fliplr(unique(masked_groups));
n_unique_groups = length(unique_groups);

% hold on;
% thresholds_std_pos = ternary_thresholds.*(ternary_stdev_ratios - 1);
% thresholds_std_neg = ternary_thresholds - ternary_thresholds./ternary_stdev_ratios;
% h = errorbar(pred_scaling*axis_scaling(mask).*predictions(mask), ...
%     axis_scaling(mask).*ternary_thresholds(mask), ...
%     axis_scaling(mask).*thresholds_std_neg(mask), ...
%     axis_scaling(mask).*thresholds_std_pos(mask), ...
%     'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
% h.CapSize = 0;

masked_predictions = predictions(mask);
masked_thresholds = ternary_thresholds(mask);
masked_axis_scaling = axis_scaling(mask);
group_preds = zeros(n_unique_groups, 1);
group_thresh = zeros(n_unique_groups, 1);
for i = 1:n_unique_groups
    sub_mask_group = strcmp(masked_groups, unique_groups{i});
    sub_preds = masked_axis_scaling(sub_mask_group).*masked_predictions(sub_mask_group);
    sub_thresh = masked_axis_scaling(sub_mask_group).*masked_thresholds(sub_mask_group);
    
    group_preds(i) = mean(sub_preds);
    group_thresh(i) = mean(sub_thresh);
end

% scatterfit(pred_scaling*group_preds, group_thresh, 'scatteropts', {'density', false}, ...
%     'fitopts', {'line', [1 0], 'style', {'k--'}});
scatterfit(pred_scaling*group_preds, group_thresh, 'scatteropts', {'density', false}, ...
    'fitopts', {'style', {'k--'}});

t_min = min(group_thresh);
t_max = max(group_thresh);

ylim([0.9*t_min 1.1*t_max]);
% axis equal;

xlabel('Predicted thresholds');
ylabel('Actual thresholds');

for i = 1:n_unique_groups
    text(pred_scaling*group_preds(i)+0.01, group_thresh(i), unique_groups{i}, ...
        'fontsize', 6);        
end

beautifygraph;
preparegraph;

%% Compare distribution in AB_1_1 and AB_1_2 directions

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 18 3];

idx1 = find(strcmp(ternary_groups, 'AB_1_1'), 1);
idx2 = find(strcmp(ternary_groups, 'AB_1_2'), 1);

% get mapping from the 99 probabilities to texture groups
mtc = processBlock('mtc', 3);
unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
groups = repelem(unique_groups, 3);

ni_ev_full_centered = ni_ev_full - 1/3;
hist_edges = linspace(-1, 1, 20);
for i0 = 1:12
    % 1 --> 7, 2 --> 6, 3 --> 5, ...
%     j0 = mod(1 - i0 + 6, 12) + 1;
    % 1 --> 1, 2 --> 12, 3 --> 11, ...
    j0 = mod(1 - i0, 12) + 1;
    i = idx1 + i0 - 1;
    j = idx2 + j0 - 1;
    
    axis1 = zeros(size(ni_ev_full_centered, 2), 1);
    mask1 = strcmp(groups, ternary_groups{i});
    axis1(mask1) = ternary_axes{i} - 1/3;
    projection1 = ni_ev_full_centered*axis1 / sum(axis1 .^ 2);
    
    axis2 = zeros(size(ni_ev_full_centered, 2), 1);
    mask2 = strcmp(groups, ternary_groups{j});
    axis2(mask2) = ternary_axes{j} - 1/3;
    projection2 = ni_ev_full_centered*axis2 / sum(axis2 .^ 2);
    
    subplot(1, 12, i0);
    [hist1, edges1] = histcounts(projection1, hist_edges);
    [hist2, edges2] = histcounts(projection2, hist_edges);
    
    bar(edges1(1:end-1), hist1);
    hold on;
    bar(edges2(1:end-1), -hist2);
    
    yl = ylim;
    
    text(0, yl(2)*2/3, [ternary_groups{idx1} '(' int2str(i - idx1 + 1) ')']);
    text(0, yl(1)*2/3, [ternary_groups{idx2} '(' int2str(j - idx2 + 1) ')']);
%     legend({['+: ' ternary_groups{idx1} '(' int2str(i - idx1 + 1) ')'], ...
%         ['-: ' ternary_groups{idx2} '(' int2str(j - idx2 + 1) ')']}, ...
%         'location', 'south', 'box', 'off', 'fontsize', 8);
    
    yrng = max(ylim);
    ylim([-yrng yrng]);
end

%% Make plots for each texture group

unique_groups = fliplr(unique(ternary_groups));

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

plotter = MatrixPlotter(length(unique_groups));
while plotter.next
    i = plotter.index;
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
    mask_idxs = find(strcmp(ternary_groups, unique_groups{i}));
    measured_points = [];
    predicted_points = [];
    h_meas = [];
    h_pred = [];
    for j = 1:length(mask_idxs)
        crt_idx = mask_idxs(j);
        crt_axis = ternary_axes{crt_idx};
        
        % measured threshold
        crt_threshold = ternary_thresholds(crt_idx);
        crt_thresh_pos = [1/3 1/3 1/3]*(1 - crt_threshold) + crt_axis*crt_threshold;
                
        % predicted threshold
        crt_predicted_threshold = pred_scaling*predictions(crt_idx);
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

%% Correlation between threshold and inverse norm of axis

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(ternary_thresholds) & ~isnan(predictions));
% exclude A_1 directions because we equalized/contrast-adapted natural
% image patches, so that direction is not properly represented by our
% analysis
% mask(strcmp(ternary_groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
% mask(isnan(stdev_ratios)) = false;

masked_groups = ternary_groups(mask);
unique_groups = fliplr(unique(masked_groups));
% put A_1 last, if it exists
if strcmp(unique_groups{1}, 'A_1')
    unique_groups = [unique_groups(2:end) unique_groups(1)];
end
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

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];
hold on;

masked_predictions = predictions(mask);
masked_thresholds = ternary_thresholds(mask);
masked_axes = ternary_axes(mask);
for i = 1:n_unique_groups
    sub_mask_group = strcmp(masked_groups, unique_groups{i});
    sub_thresh = masked_thresholds(sub_mask_group);
    sub_preds = masked_predictions(sub_mask_group);
    sub_axes = masked_axes(sub_mask_group);
    
    crt_color = 0.5 + 0.5*group_colors(i, :);
    inv_norms = cellfun(@(a) 1./norm(a)^2, sub_axes);
    smartscatter(sub_thresh, inv_norms, 'density', false, ...
        'color', crt_color);
    crt_mask = isfinite(sub_thresh) & isfinite(inv_norms);
    rho = corr(flatten(sub_thresh(crt_mask)), flatten(inv_norms(crt_mask)));
    
    crt_mask = isfinite(sub_preds) & isfinite(inv_norms);
    rho_preds = corr(flatten(sub_preds(crt_mask)), flatten(inv_norms(crt_mask)));
%     disp([unique_groups{i} ': corr(preds, inv_norms) = ' num2str(rho_preds, '%.2f') ...
%         ', p = ' num2str(pval_preds, '%.2g')]);
    
%     disp([unique_groups{i} ': corr(thresh, inv_norms) = ' num2str(rho, '%.2f') ...
%         ', corr(preds, inv_norms) = ' num2str(rho_preds, '%.2f')]);
    disp([unique_groups{i} ': corr(thresh, inv_norms) = ' num2str(rho, '%.2f')]);
end

%% Show location of tex_axes

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

for i = 1:12
    crt_axis = ternary_axes{i};
    
    % this is a litle roundabout, but easier conceptually
    crt_t1 = (3*crt_axis(2) - 1)/2;
    crt_t2 = (crt_axis(3) - crt_axis(1)) * max_t2;
    
    plot(crt_t1, crt_t2, 'ks', 'linewidth', 1)
end

% [1, 0, 0] -> [-1/2, -sqrt(3)/2]
% [0, 1, 0] -> [1, 0]
% [0, 0, 1] -> [-1/2,  sqrt(3)/2]
%
% d([1, 0, 0], [0, 1, 0]) = sqrt(3)
% d([0, 0, 1], [0, 1, 0]) = sqrt(3)
% d([1, 0, 0], [0, 0, 1]) = sqrt(3)
%
% --> equilateral triangle

% norm(t1, t2)^2 = (3*a2-1)^2/4 + 3/4*(a3 - a1)^2
%                = 3/4*(a1^2 + 3*a2^2 + a3^2) - 3/2*a1*a3 + 1/4 - 3*a2/2
%                = 3/4*norm(a)^2 + 1/4 - 3/2*(a2*(1-a2) + a1*a3)
%                = (...)
%                = 3/2*norm(a)^2 - 1/2 + 3/4*(1-sum(a))*(1 + a1 - a2 + a3)

beautifygraph;
preparegraph;

%% Get threshold predictions for a larger number of axes

n_axes_per_plane = 64;

unique_ternary_groups = fliplr(unique(ternary_groups));
many_ax_groups = repelem(unique_ternary_groups, n_axes_per_plane);

angle_vec = 0:2*pi/n_axes_per_plane:2*pi*(n_axes_per_plane-1)/n_axes_per_plane;
each_ax_set_2d = [cos(angle_vec(:)) sin(angle_vec(:))];

% t1 = (3*a2 - 1)/2
% t2 = (a3 - a1)*sqrt(3)/2
% a1 + a2 + a3 = 1
% 
% find a1, a2, a3 as a function of t1, t2
%
% a1 = (1 - t1 - sqrt(3)*t2)/3
% a2 = (1 + 2*t1)/3
% a3 = (1 - t1 + sqrt(3)*t2)/3

each_ax_set_3d = num2cell([(1 - each_ax_set_2d(:, 1) - sqrt(3)*each_ax_set_2d(:, 2))/3 ...
    (1 + 2*each_ax_set_2d(:, 1))/3 ...
    (1 - each_ax_set_2d(:, 1) + sqrt(3)*each_ax_set_2d(:, 2))/3], 2);

many_ax_axes = flatten(repmat(each_ax_set_3d, 1, length(unique_ternary_groups)));

ni_cov = cov(ni_ev_full);
% scaling_mat = sqrtm(ni_cov);
scaling_mat = sqrtm(ni_cov + 1e-10*eye(size(ni_cov)));
% scaling_mat = ni_cov;
% scaling_mat = ni_cov;

% get mapping from the 99 probabilities to texture groups
mtc = processBlock('mtc', 3);
unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
groups = repelem(unique_groups, 3);

many_ax_predictions = zeros(size(many_ax_groups));
for i = 1:length(many_ax_groups)
    crt_axis = zeros(size(scaling_mat, 1), 1);
    crt_mask = strcmp(groups, many_ax_groups{i});
    if sum(crt_mask) ~= 3
        error('Some group appears something different from 3 times.');
    end
    crt_axis(crt_mask) = many_ax_axes{i} - 1/3;
%     many_ax_predictions(i) = sum(crt_axis.^2)/(crt_axis'*scaling_mat*crt_axis);
    many_ax_predictions(i) = 1/sqrt(crt_axis(:)'*scaling_mat*crt_axis(:));
end

%% Make plots for each texture group with many axes per plane

unique_groups = fliplr(unique(many_ax_groups));

max_t2 = sqrt(3)/2;
angle_range = linspace(0, 2*pi, 100);

plotter = MatrixPlotter(length(unique_groups));
while plotter.next
    i = plotter.index;
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
    
    % draw predicted thresholds
    mask_idxs = find(strcmp(many_ax_groups, unique_groups{i}));
    predicted_points = [];
    h_pred = [];
    for j = 1:length(mask_idxs)
        crt_idx = mask_idxs(j);
        crt_axis = many_ax_axes{crt_idx};
        
        % predicted threshold
        crt_predicted_threshold = pred_scaling*many_ax_predictions(crt_idx);
        crt_pred_thresh_pos = [1/3 1/3 1/3]*(1 - crt_predicted_threshold) + ...
            crt_axis*crt_predicted_threshold;
        
        % this is a litle roundabout, but easier conceptually
        crt_pred_t1 = (3*crt_pred_thresh_pos(2) - 1)/2;
        crt_pred_t2 = (crt_pred_thresh_pos(3) - crt_pred_thresh_pos(1)) * max_t2;
        
        % store points for fitting ellipse
        if ~any(isnan([crt_pred_t1 crt_pred_t2]))
            predicted_points = [predicted_points ; [crt_pred_t1 crt_pred_t2]]; %#ok<AGROW>
        end
        
        % plot
        h_pred = [h_pred plot(crt_pred_t1, crt_pred_t2, 'r.')]; %#ok<AGROW>
    end
    
    legend(h_pred(1), 'predictions');
    
    % show ellipses
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

%% Show the thresholds in each of the planes

% % expand natural stats to contain all 99 probabilities
% ni_ev_full = reshape(natural_stats.ev, [], 2, 33);
% ni_ev_full(:, 3, :) = 1 - sum(ni_ev_full, 2);
% ni_ev_full = reshape(ni_ev_full, [], 99);
% 
% ni_ev_full_centered = ni_ev_full - 1/3;
% 
% % get mapping from the 99 probabilities to texture groups
% mtc = processBlock('mtc', 3);
% unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
% groups = repelem(unique_groups, 3);
% 
% % focus on groups that have measured thresholds
% measured_mask = isfinite(ternary_thresholds);
% measured_groups = ternary_groups(measured_mask);
% measured_axes = ternary_axes(measured_mask);
% measured_thresholds = ternary_thresholds(measured_mask);
% plotter = MatrixPlotter(length(measured_groups));
% while plotter.next
%     i = plotter.index;
%     
%     % make histogram of natural statistics
%     crt_axis_centered = zeros(size(ni_ev_full_centered, 2), 1);
%     crt_mask = strcmp(groups, measured_groups{i});
%     if sum(crt_mask) ~= 3
%         error('Some group appears something different from 3 times.');
%     end
%     crt_axis_centered(crt_mask) = measured_axes{i} - 1/3;
%     ni_projection = ni_ev_full_centered*crt_axis_centered / sum(crt_axis_centered .^ 2);
%     histogram(ni_projection, 20);
%     
%     hold all;
%     
%     % find location of threshold
%     plot([measured_thresholds(i) measured_thresholds(i)], ylim, '--r');
%     
%     title(measured_groups{i});
% end