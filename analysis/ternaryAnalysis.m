% analyze the G=3 data using the G=3 stats

%% Load ternary stats

% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog.mat'));
% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_nolog_contrastadapt_with_focus.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_contrastadapt_with_focus.mat'));
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
ni_ev_full = expandev(natural_stats.ev, 3);

%% Load the psychophysics data

ternary_avg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));
ternary_per_subject = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
    'subjects', '*', 'keepnan', false);

%% Show natural stats in PC1/2 plane

[pca_coeffs, pca_score, ~, ~, pca_explained] = pca(ni_ev_full);

figure;
smartscatter(pca_score(:, 1), pca_score(:, 2));
xlabel('PC1');
ylabel('PC2');

beautifygraph;
preparegraph;

figure;
plot(pca_explained);
xlabel('PC index');
ylabel('Varience explained (%)');

beautifygraph;
preparegraph;

%% The shape of noise

% start with patch at origin, add uncorrelated noise
% what is the covariance matrix of the resulting noise?
n_noise_samples = 10000;

% keep things reproducible
rng(5435);

noise_patch_size = 64;
base_noise_patch = quantize(rand(noise_patch_size), 3);

noise_size = 0.3;
noise_patches = arrayfun(@(p) quantize(...
    min(max(base_noise_patch + noise_size*randn(noise_patch_size), 0), 1), 3), ...
    1:n_noise_samples, 'uniform', false);

noise_ev = zeros(n_noise_samples, size(ni_ev_full, 2));
progress = TextProgress;
for i = 1:n_noise_samples
    [~, crt_ev] = processBlock(noise_patches{i}, 3);
    noise_ev(i, :) = expandev(crt_ev(:)', 3);
    if mod(i, 10) == 0
        progress.update(100*i/n_noise_samples);
    end
end
progress.done;

noise_cov = cov(noise_ev);
% normalize the noise covariance matrix so that its highest eigenvalue is 1
noise_cov = noise_cov / max(eig(noise_cov));

%% Get threshold predictions

% whether to use the noise estimate or not
use_noise = true;

% get gain matrix from efficient coding
ni_cov = cov(ni_ev_full);
if ~use_noise
    noise_cov_to_use = eye(size(ni_cov));
    in_noise_amt = 2.0;
    lag_choice = 3e-5;
else
    noise_cov_to_use = noise_cov;
    in_noise_amt = 0.5;
    lag_choice = 5e-6;
end
gain = solveLinearEfficientCoding(ni_cov, in_noise_amt*noise_cov_to_use, ...
    eye(size(ni_cov, 1)), lag_choice);

% get thresholds based on gain matrix
predictions0 = gainsToThresholds(gain, cellfun(@(v) v - 1/3, ...
    ternaryextdir(ternary_avg.groups, ternary_avg.directions), 'uniform', false));

% scale predictions to be as close as possible to measurements
[~, predictions, pred_sse] = fitscale(predictions0, ternary_avg.thresholds, ...
    'stds', diff(ternary_avg.threshold_intervals, [], 2), ...
    'exclude', strcmp(ternary_avg.groups, 'A_1'));

%% Find optimal noise size and Lagrange multiplier

tic;
% noise_sizes = linspace(0.1, 1, 50);
noise_sizes = logspace(-1, 3, 50);
lags = logspace(-10, -2, 50);
sses = zeros(length(noise_sizes), length(lags));
for i = 1:length(noise_sizes)
    crt_noise = noise_sizes(i);
    for j = 1:length(lags)
        crt_lag = lags(j);
        
        crt_gain = solveLinearEfficientCoding(ni_cov, crt_noise*noise_cov_to_use, ...
            eye(size(ni_cov, 1)), crt_lag);
        crt_preds0 = gainsToThresholds(crt_gain, cellfun(@(v) v - 1/3, ...
            ternaryextdir(ternary_avg.groups, ternary_avg.directions), 'uniform', false));

        [~, ~, sses(i, j)] = fitscale(crt_preds0, ternary_avg.thresholds, ...
            'stds', diff(ternary_avg.threshold_intervals, [], 2), ...
            'exclude', strcmp(ternary_avg.groups, 'A_1'));
    end
end
disp(['Took ' num2str(toc, '%.2f') ' seconds.']);

%% Get threshold predictions (old)

ni_cov = cov(ni_ev_full);
% scaling_mat = sqrtm(ni_cov);
scaling_mat = sqrtm(ni_cov + 1e-10*eye(size(ni_cov)));
% scaling_mat = ni_cov;

% get mapping from the 99 probabilities to texture groups
mtc = processBlock('mtc', 3);
unique_groups = cellfun(@(s) s.name, mtc.coord_groups, 'uniform', false);
groups = repelem(unique_groups, 3);

predictions_old = zeros(size(ternary_avg.groups));
for i = 1:length(ternary_avg.groups)
    crt_axis = zeros(size(scaling_mat, 1), 1);
    crt_mask = strcmp(groups, ternary_avg.groups{i});
    if sum(crt_mask) ~= 3
        error('Some group appears something different from 3 times.');
    end
    crt_axis(crt_mask) = ternary_avg.directions{i} - 1/3;
    predictions_old(i) = 1/sqrt(crt_axis'*scaling_mat*crt_axis);
end

old_mask = isfinite(ternary_avg.thresholds) & isfinite(predictions_old);
pred_scaling = dot(ternary_avg.thresholds(old_mask), predictions_old(old_mask)) / ...
    sum(predictions_old(old_mask).^2);
predictions_old = predictions_old*pred_scaling;

%% Compare predicted to actual thresholds

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];

showPredictionMatch(predictions, ternary_avg.thresholds, ternary_avg.groups, ...
    'measintervals', ternary_avg.threshold_intervals, 'exclude', strcmp(ternary_avg.groups, 'A_1'));

preparegraph;

%% Compare per-plane average predicted to actual thresholds across subjects

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 9 7];

predictions_per_subject = matchThresholds(predictions, ...
    ternary_avg.groups, ternary_avg.directions, ...
    ternary_per_subject.groups, ternary_per_subject.directions);
showPredictionMatch(predictions_per_subject, ternary_per_subject.thresholds, ...
    ternary_per_subject.groups, ...
    'subjavg', true, 'subjects', ternary_per_subject.subjects, ...
    'measintervals', ternary_per_subject.threshold_intervals, ...
    'exclude', strcmp(ternary_per_subject.groups, 'A_1'));

preparegraph;

%% Make plots for each texture group

ternaryPredictionMatchPerGroup(predictions, ternary_avg, 'ellipses', false, ...
    'exclude', strcmp(ternary_avg.groups, 'A_1'));

%% SCRATCH

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