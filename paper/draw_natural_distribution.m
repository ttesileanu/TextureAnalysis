% make plot of natural image distribution, showing in-focus region

%% Load data

natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_contrastadapt_with_focus.mat'));
% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary.mat'));
natural_stats = natural_stats0.res{1};

%% Show natural stats in PC1/2 plane

[pca_coeffs, pca_score, ~, ~, pca_explained] = pca(natural_stats.ev);

ev_mc = bsxfun(@minus, natural_stats.ev, mean(natural_stats.ev, 1));
infocus_projected = ev_mc(natural_stats.focus.clusterIds == natural_stats.focus.focusCluster, :)*...
    pca_coeffs;

figure;
hold on;
scatter(pca_score(:, 1), pca_score(:, 2), [], [0.5 0.5 0.5]);
smartscatter(infocus_projected(:, 1), infocus_projected(:, 2));
xlabel('PC1');
ylabel('PC2');

beautifygraph;
preparegraph;

%% Show natural stats in second-order -- fourth-order plane

mtc = processBlock('mtc', 3);
group_names = flatten(repmat(cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false), 2, 1));

group1 = 'AB_1_2';
group2 = 'ABCD_1_2_1_2';

% second_order_mask = cellfun(@(s) length(s) == 6, group_names);
second_order_mask = find(strcmp(group_names, group1), 1);
second_order = mean(natural_stats.ev(:, second_order_mask), 2);

% fourth_order_mask = cellfun(@(s) length(s) == 12, group_names);
fourth_order_mask = find(strcmp(group_names, group2), 1);
fourth_order = mean(natural_stats.ev(:, fourth_order_mask), 2);

focus_mask = (natural_stats.focus.clusterIds == natural_stats.focus.focusCluster);
second_order_focus = second_order(focus_mask);
fourth_order_focus = fourth_order(focus_mask);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 2.5 2.25];
hold on;
smartscatter(second_order, fourth_order, 'density', false, 'color', [0.5 0.5 0.5], ...
    'alpha', 0.05, 'maxpoints', inf, 'size', 5);
% smartscatter(second_order_focus, fourth_order_focus, 'alpha', 0.4, 'maxpoints', inf);
smartscatter(second_order_focus, fourth_order_focus, 'density', false, 'color', [0    0.4470    0.7410], ...
    'alpha', 0.05, 'maxpoints', inf, 'size', 5);
% xlabel('Second order');
% ylabel('Fourth order');
xlabel(group1);
ylabel(group2);

axis equal;

xlim([0.1 0.999]);
ylim([0.2 0.999]);

beautifygraph('fontscale', 0.6667, 'ticksize', 11);
preparegraph;

safe_print(fullfile('figs', 'draft', 'focus_locus'), 'type', 'png', 'printopts', {'-r600'});