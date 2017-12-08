% Reanalyzing binary psychophysics using continuous texture statistics.

%% Load psychophysics data, and binary image statistics

data0 = open(fullfile('data', 'PPandNIstats.mat'));

%% Load natural image stats

% natural_stats0 = open(fullfile('save', 'natural_nosky_multiscale_60_AF_1_2_3_4_5_6.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_continuous_with_focus.mat'));
natural_stats_bin0 = open(fullfile('save', 'natural_nosky_binary_with_focus.mat'));
% keep only blockAF > 1
natural_stats = natural_stats0.res(natural_stats0.N_values > 1);
natural_stats_bin = natural_stats_bin0.res(natural_stats_bin0.N_values > 1);

%% Find scalings between binary and continuous stats

% use the scalings?
use_ni_scalings = true;

ni_scalings = repmat({ones(1, 10)}, size(natural_stats));
if use_ni_scalings
    for k = 1:length(natural_stats)
        for i = 1:10
%             [~, fit_stats] = drawfitline(...
%                 natural_stats{k}.ev(:, i), ...
%                 natural_stats_bin{k}.ev(:, i), ...
%                 'nodraw', true, 'intercept', 0);
%             ni_scalings{k}(i) = fit_stats.a;
            ni_scalings{k}(i) = std(natural_stats_bin{k}.ev(:, i)) / ...
                std(natural_stats{k}.ev(:, i));
        end
    end
end

%% Scale the continuous stats to look more like the binary stats

natural_stats_scaled = natural_stats;
for k = 1:length(natural_stats_scaled)
    natural_stats_scaled{k}.ev = bsxfun(@times, ...
        natural_stats{k}.ev, ni_scalings{k});
    natural_stats_scaled{k}.covM = cov(natural_stats_scaled{k}.ev);
end

%% How do my binary statistics compare to Ann's?

ann_ev = data0.dataNI.indA(6).ev;
% my_ev = natural_stats_bin{1}.ev;
my_ev = natural_stats_bin{1}.ev(natural_stats_bin{1}.focus.clusterIds == ...
    natural_stats_bin{1}.focus.focusCluster, :);

c_labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

plot_pairs = [1 2 ; 3 4 ; 5 6 ; 7 10 ; 1 10 ; 4 10 ; 1 7 ; 6 9];

fig = figure;
fig.Position = [fig.Position([1 2]) 3*fig.Position(3) fig.Position(4)];
for k = 1:size(plot_pairs, 1)
    subplot(2, 4, k);
    
    hold on;
    
    i = plot_pairs(k, 1);
    j = plot_pairs(k, 2);
    
    smartscatter(ann_ev(:, i), ann_ev(:, j), 'color', [1 0 0], 'density', false, 'alpha', 0.1);
    smartscatter(my_ev(:, i), my_ev(:, j), 'color', [0, 0, 1], 'density', false, 'alpha', 0.1);
    
    xlabel(c_labels{i});
    ylabel(c_labels{j});
    
    legend({'new', 'elife'});
    
    beautifygraph;
end
preparegraph;

%% Predict thresholds

% ni_covs = cellfun(@(s) s.covM, natural_stats_scaled, 'uniform', false);
ni_covs = cellfun(@(s) cov(s.ev(s.focus.clusterIds == s.focus.focusCluster, :)), ...
    natural_stats_scaled, 'uniform', false);
ni_precision = cellfun(@(natural_cov) inv(natural_cov), ni_covs, 'uniform', false);

ni_sensitivities = cellfun(@(m) sqrt(1000*diag(m)), ni_covs, 'uniform', false);
ni_thresholds = cellfun(@(t) 1./t, ni_sensitivities, 'uniform', false);

% collect predictions from binary statistics

% focus on only some of the binary results, as in Hermunstad et al.
ann_mask = arrayfun(@(s) ismember(s.N, [2, 4]) && ismember(s.R, [32, 48, 64]), ...
    data0.dataNI.indA);
% ann_covs = arrayfun(@(s) s.covM, data0.dataNI.indA(ann_mask), 'uniform', false);
ann_covs = arrayfun(@(s) cov(s.ev), data0.dataNI.indA(ann_mask), 'uniform', false);
ann_precision = cellfun(@(natural_cov) inv(natural_cov), ann_covs, 'uniform', false);

ann_sensitivities = cellfun(@(m) sqrt(1000*diag(m)), ann_covs, 'uniform', false);
ann_thresholds = cellfun(@(t) 1./t, ann_sensitivities, 'uniform', false);

binary_covs = cellfun(@(s) cov(s.ev(s.focus.clusterIds == s.focus.focusCluster, :)), ...
    natural_stats_bin, 'uniform', false);
binary_precision = cellfun(@(natural_cov) inv(natural_cov), binary_covs, 'uniform', false);

binary_sensitivities = cellfun(@(m) sqrt(1000*diag(m)), binary_covs, 'uniform', false);
binary_thresholds = cellfun(@(t) 1./t, binary_sensitivities, 'uniform', false);

%% Transform psychophysics data to thresholds and sensitivities

pp_covs = arrayfun(@(s) s.covM, data0.dataPP.indS, 'uniform', false);
pp_precision = cellfun(@(natural_cov) inv(natural_cov), pp_covs, 'uniform', false);

% the 1000 factor was used for plotting in Hermunstad et al.
pp_sensitivities = arrayfun(@(s) sqrt(1000*diag(s.covM)), data0.dataPP.indS, 'uniform', false);
pp_thresholds = cellfun(@(t) 1./t, pp_sensitivities, 'uniform', false);

%% Make threshold plot (Fig. 3A in Hermunstad et al. 2014)

% find average pp sensitivities
pp_avg_diagsensitivity = mean(cell2mat(pp_sensitivities), 2);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 12 8];
hold on;

% draw the predictions from natural statistics

c_labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};
% locations of the 9 statistics (matching c_labels except gamma) on the x axis
all_x_positions = [1:2 5:6 9:12 15];

all_ni_diag = [ni_sensitivities(:)' binary_sensitivities(:)' ann_sensitivities(:)'];
all_ni_shifts = [zeros(1, length(ni_sensitivities)) 0.2*ones(1, length(binary_sensitivities)) ...
    0.4*ones(1, length(ann_sensitivities))];

select_half = @(m) m(1:end/2, :);
all_ni_colors = [...
    cool(length(ni_sensitivities)) ; ...
    cool(length(binary_sensitivities)) ; ...
    cool(length(ann_sensitivities))];

continuous_ni_labels = cellfun(@(s) ['continuous N=' int2str(s.options.blockAF) ...
    ', R=' int2str(s.options.patchSize(1))], natural_stats_scaled, 'uniform', false);
ann_ni_labels = arrayfun(@(s) ['eLife N=' int2str(s.N) ', R=' int2str(s.R)], ...
    data0.dataNI.indA(ann_mask), 'uniform', false);
binary_ni_labels = cellfun(@(s) ['binary N=' int2str(s.options.blockAF) ...
    ', R=' int2str(s.options.patchSize(1))], natural_stats_bin, 'uniform', false);
all_ni_labels = [continuous_ni_labels(:)' binary_ni_labels(:)' ann_ni_labels(:)'];

symbol_size = 7;

for i = 1:length(all_ni_diag)
    crt_sensitivities = all_ni_diag{i};
    
    % find optimal scaling
    crt_scaling = dot(crt_sensitivities(2:end), pp_avg_diagsensitivity(2:end)) / ...
        norm(crt_sensitivities(2:end))^2;
    
    plot(all_x_positions + all_ni_shifts(i), crt_scaling*crt_sensitivities(2:end), 'o', ...
        'markersize', symbol_size, 'color', all_ni_colors(i, :), 'linewidth', 2);
end

% draw the psychophysics results
pp_colors = select_half(hot(2*length(pp_sensitivities)));
pp_labels = arrayfun(@(s) ['data ' s.subjID], data0.dataPP.indS, 'uniform', false);
measurement_x_positions = [3, 7, 13, 16];
measurement_subselect_idxs = [2, 4, 6, 10];
for i = 1:length(pp_sensitivities)
    crt_sensitivities = pp_sensitivities{i};
    plot(measurement_x_positions, crt_sensitivities(measurement_subselect_idxs), ...
        's', 'color', pp_colors(i, :), 'markersize', symbol_size, 'linewidth', 2);
end

% set up labels and ticks
title('Texture sensitivity from natural images vs. psychophysics');
ylabel('Sensitivity');

xlim([0 17]);
ylim([0 4.5]);

legend([all_ni_labels pp_labels]);

beautifygraph;

set(gca, 'xtick', all_x_positions, 'xticklabel', c_labels(2:end), 'xminortick', 'off');

preparegraph;

print('-dpdf', fullfile('figs', 'sensitivity_continuous_binary_experiment.pdf'));

%% Make ellipse plot (Fig. 3B in Hermunstad et al. 2014)

ni_covs_scaled = cell(size(ni_covs));

for i = 1:length(ni_covs)
    crt_sensitivities = ni_sensitivities{i};
    
    % find optimal scaling
    crt_scaling = dot(crt_sensitivities(2:end), pp_avg_diagsensitivity(2:end)) / ...
        norm(crt_sensitivities(2:end))^2;
    
    ni_covs_scaled{i} = ni_covs{i}*crt_scaling^2;
end

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 12 12];
hold on;

binary_chosen_ellipse = 1;
ni_chosen_ellipse = 1;
ann_chosen_ellipse = 1;
% binary_mask_idxs = find(binary_mask);
all_ellipse_covs = [...
    pp_covs ...
    binary_covs(binary_chosen_ellipse) ...
    ann_covs(ann_chosen_ellipse) ...
    ni_covs_scaled(ni_chosen_ellipse) ...
];
% all_ellipse_precision = {};
all_ellipse_colors = [pp_colors ; 0 0.8 0 ; 0 0.2 0.3 ; 0 0 1];
all_ellipse_line_width = [0.5*ones(1, length(pp_precision)) 1 1 1];
all_ellipse_labels = [...
    pp_labels ...
    binary_ni_labels(binary_chosen_ellipse) ...
    ann_ni_labels(ann_chosen_ellipse) ...
    arrayfun(@(i) ...
        ['continuous N=' int2str(natural_stats_scaled{i}.options.blockAF) ...
        ' R=' int2str(natural_stats_scaled{i}.options.patchSize(1))], ni_chosen_ellipse(:)', 'uniform', false) ...
];

for i = 2:10
    for j = 2:i-1
        % figure out scaling (such that the maximum major axis is set to 1)
        all_a2 = cellfun(@(m) max(1./eig(m([j i], [j i]))), all_ellipse_covs);
        cov_scaling = 0.25/max(all_a2);
        
        % draw ellipses
        for k = 1:length(all_ellipse_covs)
            crt_cov = all_ellipse_covs{k};
            crt_mat = inv(crt_cov([j i], [j i]));
            
            ellipse(j-1, 10-i+1, cov_scaling*crt_mat, ...
                'color', all_ellipse_colors(k, :), ...
                'linewidth', all_ellipse_line_width(k)); %#ok<MINV>
        end
    end
end

axis equal;
xlim([0 10]);
ylim([0 10]);

set(gca, 'xtick', (1:9), 'xticklabel', c_labels(2:end));
set(gca, 'ytick', (1:9), 'yticklabel', fliplr(c_labels(2:end)));
% set(gca, 'ydir', 'reverse');

title('Precision matrix comparisons between natural images and psychophysics');

legend(all_ellipse_labels);

beautifygraph;

set(gca, 'xminortick', 'off', 'yminortick', 'off');

preparegraph;

print('-dpdf', fullfile('figs', 'ellipses_continuous_binary_experiment.pdf'));