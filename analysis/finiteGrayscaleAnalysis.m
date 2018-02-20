% analyze the G=3 data using the G=3 stats

%% Load ternary stats

% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog_contrastadapt.mat'));
natural_stats = natural_stats0.res{1};

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
    predictions(i) = sum(crt_axis.^2)/(crt_axis'*scaling_mat*crt_axis);
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
pred_scaling = dot(ternary_thresholds(mask), predictions(mask)) / ...
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
h = errorbar(pred_scaling*predictions(mask), ternary_thresholds(mask), ...
    thresholds_std_neg(mask), thresholds_std_pos(mask), ...
    'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
h.CapSize = 0;

masked_predictions = predictions(mask);
masked_thresholds = ternary_thresholds(mask);
for i = 1:n_unique_groups
    sub_mask_group = strcmp(masked_groups, unique_groups{i});
    sub_preds = masked_predictions(sub_mask_group);
    sub_thresh = masked_thresholds(sub_mask_group);
    
    crt_color = 0.5 + 0.5*group_colors(i, :);
    
    for j = 1:length(sub_preds)
        for k = j+1:length(sub_preds)
            line(pred_scaling*sub_preds([j k]), sub_thresh([j k]), 'linewidth', 0.5, ...
                'color', crt_color);
        end
    end
end

th_h = smartscatter(pred_scaling*predictions(mask), ternary_thresholds(mask), 'color', threshold_colors, ...
    'density', false);

t_min = min(ternary_thresholds(mask));
t_max = max(ternary_thresholds(mask));

ylim([0.8*t_min 1.2*t_max]);
axis equal;

% drawfitline(predicted_thresholds(mask), thresholds(mask), 'showci', false, ...
%     'corrtype', 'spearman');
drawfitline(pred_scaling*predictions(mask), ternary_thresholds(mask), 'line', [1 0], ...
    'style', {'--k'}, 'legendloc', 'northwest');
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

for i = 1:length(ternary_thresholds)
    if ~mask(i)
        continue;
    end
    
    text(pred_scaling*predictions(i)+0.01, ternary_thresholds(i), ternary_groups{i}, ...
        'fontsize', 6);        
%         [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']'], ...
end

beautifygraph;
preparegraph;

%% Compare distribution in A_B_1_1 and AB_1_2 directions

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
    j0 = mod(1 - i0 + 6, 12) + 1;
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
    
    disp([unique_groups{i} ': corr(thresh, inv_norms) = ' num2str(rho, '%.2f') ...
        ', corr(preds, inv_norms) = ' num2str(rho_preds, '%.2f')]);
end
    
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