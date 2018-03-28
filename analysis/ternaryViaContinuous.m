% Analyzing three grayscale level-data using continuous texture statistics.

%% Load psychophysics data

ternary_avg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));
ternary_per_subject = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
    'subjects', '*', 'keepnan', false);

%% Generate patches in the directions used in the psychophysics

% keep this reproducible
rng(34235);

[mapping, mapping_details] = generateTextureMapping(...
    @(patch) getnth(2, @processBlock, patch, inf), ternary_avg.groups, ternary_avg.directions);

%% Save results

save(fullfile('save', 'ternary_via_continuous.mat'), ...
    'ternary_avg', 'ternary_per_subject', 'mapping', 'mapping_details');

%% Load results

load(fullfile('save', 'ternary_via_continuous.mat'));

%% Load natural image stats

natural_stats0 = open(fullfile('save', 'natural_nosky_continuous_with_focus_contrastadapt.mat'));
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

save_tag = [save_tag sprintf('N%d_R%d', natural_stats.options.blockAF, ...
    natural_stats.options.patchSize(1))];

% fit a Gaussian to the natural image data
natural_mu = mean(natural_stats.ev, 1);
natural_cov = cov(natural_stats.ev);

%% Map ternary thresholds to continuous texture space

pp_threshold_locations = arrayfun(@(i) mapping{i}.function(ternary_avg.thresholds(i)), ...
    1:length(mapping), 'uniform', false);

%% Plot the thresholds in some 10d coordinate planes

% plot the natural image distribution in the 10d continuous texture space,
% with the mapped ternary thresholds shown on top
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

pairs = {[2 10], [2 3], [3 5], [6 8], [6 10], [4 7]};

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 12 6];

threshold_locations_mat = cell2mat(pp_threshold_locations');

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

% safe_print(fullfile('figs', 'G3vsCont', ['thresholds_in_10d' save_tag]), 'png');

%% Fit natural image predictions to psychophysics measurements

% for now just tuning the threshold radius in transformed coordinates
% (i.e., after gain predicted by efficient coding)

% start by finding a good gain matrix
in_noise_amt = 1;
lag_choice = 1e-6;

% for now use uniform input noise
noise_cov_to_use = eye(size(natural_cov));

% get optimal gain matrix based on natural image statistics
gain = solveLinearEfficientCoding(natural_cov, in_noise_amt*noise_cov_to_use, ...
    eye(size(natural_cov, 1)), lag_choice);

% given the gain and a test radius^2, what are the threshold predictions?
% [predictions, predicted_locations] = mapIntersectEllipsoid(mapping, gain'*gain, 0.1);
[predictions, prediction_locations] = optimizePredictedThresholds(...
    ternary_avg.thresholds, mapping, gain'*gain, [0.1 5], 'stds', ternary_avg.threshold_intervals, ...
    'exclude', strcmp(ternary_avg.groups, 'A_1'));

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