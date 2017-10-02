% Analyzing three grayscale level-data using continuous texture statistics.

%% Load psychophysics data

f = fopen('three_grayscale_thresholds.csv');
data0 = textscan(f, '%s', 'delimiter', '');
fclose(f);

% process the data a little

% split at delimiters
groups = strsplit(data0{1}{1}, ',', 'CollapseDelimiters', false);
% copy group name over to empty fields
last_type = '';
for i = 1:length(groups)
    if isempty(groups{i})
        groups{i} = last_type;
    else
        % update format to match that required by the texture generation
        % routines
        digit_mask = find(ismember(groups{i}, '0':'9'));
        if ~isempty(digit_mask)
            letter_part = groups{i}(1:digit_mask(1)-1);
            number_part = cell2mat(arrayfun(...
                @(c) ['_' c], groups{i}(digit_mask(1):end), 'uniform', false));
            groups{i} = [letter_part number_part];
        end
        last_type = groups{i};
    end
end
tex_axes0 = strsplit(data0{1}{2}, ',', 'CollapseDelimiters', false);
% copy axes names over to empty fields
last_axis = '';
for i = 1:length(tex_axes0)
    if isempty(tex_axes0{i})
        tex_axes0{i} = last_axis;
    else
        last_axis = tex_axes0{i};
    end
end
tex_axes = cellfun(@(s) arrayfun(@str2double, s(2:end-1)), tex_axes0, 'uniform', false);
directions = strsplit(data0{1}{3}, ',', 'CollapseDelimiters', false);
thresholds = strsplit(data0{1}{4}, ',', 'CollapseDelimiters', false);
stdev_ratios = strsplit(data0{1}{5}, ',', 'CollapseDelimiters', false);
if length(groups) ~= length(tex_axes) || length(groups) ~= length(thresholds) || ...
        length(groups) ~= length(directions) || length(groups) ~= length(stdev_ratios)
    error('Number of elements is not consistent between rows!');
end

% get rid of row names
groups = groups(2:end);
tex_axes = tex_axes(2:end);
directions = directions(2:end);
thresholds = cellfun(@str2double, thresholds(2:end));
stdev_ratios = cellfun(@str2double, stdev_ratios(2:end));

%% Generate patches in the directions used in the psychophysics

% number of patches to generate for each set of coordinates
n_per_dot = 24;

% locations of dots to generate in each direction
dot_locs = linspace(-1/2, 1, 16);

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
% each direction appears twice: once for positive, and once for negative,
% deviation from center; so we skip every other one
for i = 1:2:length(tex_axes)
    percentage = (i-1)/length(tex_axes)*100;
    disp(['Working on ' groups{i} ' ' tex_axes0{i} '... (' int2str(round(percentage)) '% done, ' ...
        num2str(toc(t0), '%.2f') ' seconds elapsed)']);
    
    crt_group = groups{i};
    crt_axis = tex_axes{i};
   
    % look at locations dot_locs in this direction
    three_g_to_cont_stats{i} = cell(1, length(dot_locs));
    for j = 1:length(dot_locs)
        crt_t = dot_locs(j);
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
    three_g_to_cont_stats{i+1} = three_g_to_cont_stats{i};
end
disp('Done.');

%% Save results

save(fullfile('save', 'three_psycho_via_continuous.mat'), ...
    'directions', 'dot_locs', 'groups', 'n_per_dot', 'patch_size', ...
    'tex_axes', 'three_g_to_cont_stats', 'thresholds', 'stdev_ratios');

%% Load results

load(fullfile('save', 'three_psycho_via_continuous.mat'));

%% Load natural image stats

natural_stats0 = open(fullfile('save', 'natural_nosky_multiscale_60_AF_1_2_3_4_5_6.mat'));
blockAF_choice = 2;
natural_stats = natural_stats0.res{blockAF_choice};

%% Map G=3 thresholds to continuous space

% fit a Gaussian to the data
natural_mu = mean(natural_stats.ev, 1);
natural_cov = cov(natural_stats.ev);
natural_sqrtcov = sqrtm(natural_cov);

threshold_locations = cell(1, length(thresholds));
ellipsoid_sizes = zeros(1, length(thresholds));

% use linear interpolation or extrapolation
for i = 1:length(thresholds)
    % thresholds can point in either direction
    crt_thresh = (-1)^(strcmp(directions{i}, 'neg'))*thresholds(i);
    
    if isnan(crt_thresh)
        threshold_locations{i} = [];
        ellipsoid_sizes(i) = nan;
        continue;
    end
    
    % can we interpolate, or do we have to extrapolate?
    k = find(dot_locs(1:end-1) <= crt_thresh & dot_locs(2:end) > crt_thresh, 1);
    if ~isempty(k)
        % yes, we can interpolate
        loc1 = mean(three_g_to_cont_stats{i}{k}, 1);
        loc2 = mean(three_g_to_cont_stats{i}{k+1}, 1);
        
        alpha = (crt_thresh - dot_locs(k)) / (dot_locs(k+1) - dot_locs(k));
        thresh_loc = (1 - alpha)*loc1 + alpha*loc2;
        
        threshold_locations{i} = thresh_loc;
    else
        % we need to extrapolate
        if strcmp(directions{i}, 'neg')
            loc1 = mean(three_g_to_cont_stats{i}{1}, 1);
            loc2 = mean(three_g_to_cont_stats{i}{2}, 1);
            der = (loc2 - loc1) / (dot_locs(2) - dot_locs(1));
            thresh_loc = loc1 + der*(crt_thresh - dot_locs(1));
            threshold_locations{i} = thresh_loc;
        else
            loc1 = mean(three_g_to_cont_stats{i}{end-1}, 1);
            loc2 = mean(three_g_to_cont_stats{i}{end}, 1);
            der = (loc2 - loc1) / (dot_locs(end) - dot_locs(end - 1));
            thresh_loc = loc2 + der*(crt_thresh - dot_locs(end));
            threshold_locations{i} = thresh_loc;
        end
    end
    
    ellipsoid_sizes(i) = threshold_locations{i}(:)'*natural_sqrtcov*threshold_locations{i}(:);
end

%% Threshold in transformed continuous space --> thresholds in G=3 space

% this is actually radius *squared*
ellipsoid_size = mean(ellipsoid_sizes, 'omitnan');

%[predicted_thresholds, predicted_thresh_locs] = predictThreeGThresholds(...
%    ellipsoid_size, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs);

%% Optimize match between predicted and measured thresholds

ignore_mask = true(size(thresholds));
ignore_mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
ignore_mask(isnan(stdev_ratios)) = false;

optim_options = optimset('display', 'iter', 'tolx', 1e-10);
normdiff_mask = @(a, b, mask) norm((a(mask) - b(mask)) / (sum(mask) - 1));
normdiff_skipinf = @(a, b) normdiff_mask(a, b, ignore_mask & isfinite(a) & isfinite(b));
spcorr_mask = @(a, b, mask) corr(a(mask), b(mask), 'type', 'spearman');
spcorr_skipinf = @(a, b) spcorr_mask(a(:), b(:), ignore_mask & isfinite(a) & isfinite(b));
[optimal_size, optimal_norm, ~, optim_details] = fminbnd(@(sz) normdiff_skipinf(...
    predictThreeGThresholds(sz, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs), ...
    thresholds), 0, 2*max(ellipsoid_sizes), optim_options);
% [optimal_size, optimal_norm, ~, optim_details] = fminbnd(@(sz) -spcorr_skipinf(...
%     predictThreeGThresholds(sz, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs), ...
%     thresholds), 0, 2*max(ellipsoid_sizes), optim_options);
[predicted_thresholds, predicted_thresh_locs] = predictThreeGThresholds(...
    optimal_size, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs);

%% Make histograms of natural image stats along rays

sigma = 0.1;
natural_stats_hist = cell(size(three_g_to_cont_stats));

n_plots = length(three_g_to_cont_stats)/2;
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

for i = 1:2:length(tex_axes)
    natural_stats_hist{i} = zeros(1, length(dot_locs));
    for j = 1:length(dot_locs)
        mid_point = mean(three_g_to_cont_stats{i}{j}, 1);
        distances = bsxfun(@minus, natural_stats.ev, mid_point);
        counts = exp(-0.5*sum((distances/sigma) .^ 2, 2));
        natural_stats_hist{i}(j) = sum(counts);
    end

    ax = axes;
    ax.Units = 'normalized';
    i_crt = (i-1)/2 + 1;
    ax.OuterPosition = [mod(i_crt - 1, nx)/nx ...
        floor((i_crt - 1) / nx)/ny 1/nx 1/ny];
    
    hold on;
    fill([dot_locs(:) ; flipud(dot_locs(:))], ...
        [natural_stats_hist{i}(:) ; zeros(length(dot_locs), 1)], [0.7 0.7 0.7]);
    
    lims = ylim;
    for k = 0:1
        crt_thresh = (-1)^k*thresholds(i + k);
        crt_std_ratio = stdev_ratios(i + k);
        
        crt_t1 = crt_thresh / crt_std_ratio;
        crt_t2 = crt_thresh * crt_std_ratio;
        
        h_patch = fill([crt_t1 crt_t1 crt_t2 crt_t2], [lims(1) lims(2) lims(2) lims(1)], ...
            [0.7 0.2 0.2]);
        h_patch.FaceAlpha = 0.12;
        
        plot([crt_thresh crt_thresh], [lims(1) lims(2)], 'r--');
    end
    
    xlim([-2, 2]);
end

%% Compare predicted to actual thresholds

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(thresholds) & ~isnan(predicted_thresholds));
% exclude A_1 directions because we equalized natural image patches, so
% that direction is not properly represented by our analysis
mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
mask(isnan(stdev_ratios)) = false;

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
smartscatter(predicted_thresholds(mask), thresholds(mask), 200, [1 0 0], '.');
drawfitline(predicted_thresholds(mask), thresholds(mask), 'showci', false, ...
    'corrtype', 'spearman');
drawfitline(predicted_thresholds(mask), thresholds(mask), 'line', [1 0], ...
    'style', {'--k'});
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

dirmap = containers.Map({'pos', 'neg'}, {'+', '-'});
for i = 1:length(thresholds)
    if ~mask(i)
        continue;
    end
    text(predicted_thresholds(i)+0.01, thresholds(i), ...
        [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']' ...
        dirmap(directions{i})], ...
        'fontsize', 6);
end

beautifygraph;
preparegraph;

print('-dpdf', fullfile('figs', ['G3_to_cont_pred_vs_meas_' int2str(blockAF_choice) ...
    '.pdf']));

%% Compare predicted to actual thresholds, separate fits by order

% exclude directions for which the psychophysics yielded no thresholds
mask = (~isnan(thresholds) & ~isnan(predicted_thresholds));
% exclude A_1 directions because we equalized natural image patches, so
% that direction is not properly represented by our analysis
mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
mask(isnan(stdev_ratios)) = false;

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
smartscatter(predicted_thresholds(mask), thresholds(mask), 200, [1 0 0], '.');
drawfitline(predicted_thresholds(low_order_mask), thresholds(low_order_mask), ...
    'showci', false, 'corrtype', 'spearman', 'legendloc', 'northwest');
drawfitline(predicted_thresholds(high_order_mask), thresholds(high_order_mask), ...
    'showci', false, 'corrtype', 'spearman', 'legendloc', 'northeast');
drawfitline(predicted_thresholds(mask), thresholds(mask), 'line', [1 0], ...
    'style', {'--k'});
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

dirmap = containers.Map({'pos', 'neg'}, {'+', '-'});
for i = 1:length(thresholds)
    if ~mask(i)
        continue;
    end
    text(predicted_thresholds(i)+0.01, thresholds(i), ...
        [groups{i} '[' arrayfun(@int2str, tex_axes{i}) ']' ...
        dirmap(directions{i})], ...
        'fontsize', 6);
end

beautifygraph;
preparegraph;

%print('-dpdf', fullfile('figs', 'G3_to_cont_pred_vs_meas_two_fits.pdf'));

%% Make some plots

% if true, instead of saving some of the scatterplots to file, all of them
% are shown, while waiting for a keypress after each plot
show_all = false;

if ~show_all
    to_draw_stats = [1, 3, 7, 9, 31, 33, 49];
else
    to_draw_stats = 1:2:length(three_g_to_cont_stats);
end

for i = 1:length(to_draw_stats)
    crt_i = to_draw_stats(i);
    
    scatterTextureStats({three_g_to_cont_stats{crt_i}, dot_locs}, '.', 'alpha', 0.2);
    
    ax = axes;
    ax.Units = 'normalized';
    ax.Position = [0 0.85 1 0.15];
    axis off;
    
    text(0.5, 0.5, ['Sweep from t=-1/2 (blue) to t=1 (red) in direction ' groups{crt_i} ...
        '[' arrayfun(@int2str, tex_axes{crt_i}) ']'], ...
        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
        'fontsize', 12, 'fontweight', 'bold');
    
    crt_group = groups{crt_i}(groups{crt_i} ~= '_');
    crt_ax = arrayfun(@int2str, tex_axes{crt_i});
    filename = ['sweep_G3_in_contspace_' crt_group '_' crt_ax '.pdf'];
    if ~show_all
        print('-dpdf', fullfile('figs', filename));
    else     
        waitforbuttonpress;
        close;
    end
end

%% Show just the means for each location

% if true, instead of saving some of the scatterplots to file, all of them
% are shown, while waiting for a keypress after each plot
show_all = true;

if ~show_all
    to_draw_stats = [1, 3, 7, 9, 31, 33, 49];
else
    to_draw_stats = 1:2:length(three_g_to_cont_stats);
end

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

pairs = {[1 10], [2 3], [4 5], [6, 7], [8 9]};

colors0 = redblue(255);
colors0 = max(min(colors0 - (sum(colors0, 2) - 1)/2, 1), 0);
cidxs = 1 + round(0.5*(1 + dot_locs)*(size(colors0, 1) - 1));
colors = colors0(cidxs, :);

for i = 1:length(to_draw_stats)
    crt_i = to_draw_stats(i);
    crt_stats = three_g_to_cont_stats{crt_i};
    crt_means = cell2mat(cellfun(@(m) mean(m, 1), crt_stats', 'uniform', false));
    crt_stds = cell2mat(cellfun(@(m) std(m, [], 1), crt_stats', 'uniform', false));
    
    fig = figure;
    fig.Units = 'inches';
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

        smartscatter(crt_means(:, idx1), crt_means(:, idx2), [], colors, '.');
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
    
    text(0.5, 0.5, ['Sweep from t=-1/2 (blue) to t=1 (red) in direction ' groups{crt_i} ...
        '[' arrayfun(@int2str, tex_axes{crt_i}) ']'], ...
        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
        'fontsize', 12, 'fontweight', 'bold');
    
    crt_group = groups{crt_i}(groups{crt_i} ~= '_');
    crt_ax = arrayfun(@int2str, tex_axes{crt_i});
%    filename = ['sweep_G3_in_contspace_' crt_group '_' crt_ax '.pdf'];
    if ~show_all
%        print('-dpdf', fullfile('figs', filename));
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
    for j = 1:length(mask_idxs)
        crt_idx = mask_idxs(j);
        crt_axis = tex_axes{crt_idx};
        crt_direction = directions{crt_idx};
        
        % measured threshold
        crt_threshold = (-1)^strcmp(crt_direction, 'neg')*thresholds(crt_idx);
        crt_thresh_pos = [1/3 1/3 1/3]*(1 - crt_threshold) + crt_axis*crt_threshold;
                
        % predicted threshold
        crt_predicted_threshold = (-1)^strcmp(crt_direction, 'neg')*predicted_thresholds(crt_idx);
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
        plot(crt_t1, crt_t2, 'kx');
        plot(crt_pred_t1, crt_pred_t2, 'r.');
    end
    
    % show ellipses
    if size(measured_points, 1) > 2
        measured_M = fit_ellipse(measured_points);
        [crt_V, crt_D] = eig(measured_M);
        crt_y = [1/sqrt(crt_D(1, 1))*cos(angle_range(:)) ...
            1/sqrt(crt_D(2, 2))*sin(angle_range(:))]';
        crt_x = crt_V*crt_y;
        plot(crt_x(1, :), crt_x(2, :), 'color', [0.5 0.5 0.5]);
    end
    if size(predicted_points, 1) > 2
        predicted_M = fit_ellipse(predicted_points);
        [crt_V, crt_D] = eig(predicted_M);
        if all(diag(crt_D) > -eps)
            crt_y = [1/sqrt(crt_D(1, 1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2, 2))*sin(angle_range(:))]';
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

print('-dpdf', fullfile('figs', 'G3_to_contspace_ellipses.pdf'));

%% Check how much thresholds change with noise radius

sizes = linspace(ellipsoid_size/50, 3*ellipsoid_size, 100);

thresholds_vs_size = arrayfun(@(sz) ...
    predictThreeGThresholds(sz, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs), ...
    sizes, 'uniform', false);

% some useful definitions
ignore_mask = true(size(thresholds));
ignore_mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
ignore_mask(isnan(stdev_ratios)) = false;

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
ylim([0.6 0.9]);
ylabel('Correlation to experiment');

yyaxis right;
plot(sizes, norm_vs_size, 'r');
ylabel('RMS from experiment');

legend({'Rank corr', 'Linear corr', 'RMS'}, 'location' ,'southeast');

xlabel('Noise radius r');
title('All directions included');

beautifygraph;
preparegraph;

print('-dpdf', fullfile('figs', 'G3_to_contspace_noiserad_dependence.pdf'));

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
ylim([0.6 0.9]);
plot(sizes, spcorr_vs_size_low_order, 'b');
plot(sizes, corr_vs_size_low_order, 'b--');
ylabel('Correlation to experiment');

yyaxis right;
plot(sizes, norm_vs_size_low_order, 'r');
ylabel('RMS from experiment');

legend({'Rank corr', 'Linear corr', 'RMS'}, 'location' ,'southeast');

xlabel('Noise radius r');
title('Focusing on second-order correlations');

beautifygraph;
preparegraph;

print('-dpdf', fullfile('figs', 'G3_to_contspace_noiserad_dependence_low_order.pdf'));

%% Get optimal based on low-order

ignore_mask = true(size(thresholds));
ignore_mask(strcmp(groups, 'A_1')) = false;
% exclude directions for which the psycophysics has infinite error bars
ignore_mask(isnan(stdev_ratios)) = false;

optim_options = optimset('display', 'iter', 'tolx', 1e-10);
normdiff_low_order = @(pred, threholds) normdiff_mask(pred, thresholds, ...
    low_order_mask & isfinite(pred) & isfinite(thresholds));
optimal_size_low_order = fminbnd(@(sz) normdiff_low_order(...
    predictThreeGThresholds(sz, three_g_to_cont_stats, directions, natural_sqrtcov, dot_locs), ...
    thresholds), 0, 3*max(ellipsoid_sizes), optim_options);
[predicted_thresholds_low_order, predicted_thresh_locs_low_order] = ...
    predictThreeGThresholds(optimal_size_low_order, three_g_to_cont_stats, ...
    directions, natural_sqrtcov, dot_locs);

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
        plot(crt_t1, crt_t2, 'kx');
        plot(crt_pred_t1, crt_pred_t2, 'r.');
    end
    
    % show ellipses
    if size(measured_points, 1) > 2
        measured_M = fit_ellipse(measured_points);
        [crt_V, crt_D] = eig(measured_M);
        crt_y = [1/sqrt(crt_D(1, 1))*cos(angle_range(:)) ...
            1/sqrt(crt_D(2, 2))*sin(angle_range(:))]';
        crt_x = crt_V*crt_y;
        plot(crt_x(1, :), crt_x(2, :), 'color', [0.5 0.5 0.5]);
    end
    if size(predicted_points, 1) > 2
        predicted_M = fit_ellipse(predicted_points);
        [crt_V, crt_D] = eig(predicted_M);
        if all(diag(crt_D) > -eps)
            crt_y = [1/sqrt(crt_D(1, 1))*cos(angle_range(:)) ...
                1/sqrt(crt_D(2, 2))*sin(angle_range(:))]';
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

print('-dpdf', fullfile('figs', 'G3_to_contspace_ellipses_low_order.pdf'));