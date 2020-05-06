% Find directions with strongest mismatches between NI predictions and
% psychophysical measurements.

%% Load some colors

[colors, colorDict] = get_palette();

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the NI predictions

load(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32_square.mat'));

%% Sanity checks

assert(isequal(ternaryAvg.groups, predictions.groups));
assert(isequal(ternaryAvg.directions, predictions.directions));

%% Focus on relevant data

% we know we'll not get the first-order stats right, so get rid of those
group_mask = ~strcmp(ternaryAvg.groups, 'A_1');

%% Find absolute and log errors, and log errors divided by uncertainties

absErrors = predictions.thresholds(group_mask) - ternaryAvg.thresholds(group_mask);
logErrors = log(predictions.thresholds(group_mask)) - ...
    log(ternaryAvg.thresholds(group_mask));
logSigmas = diff(log(ternaryAvg.thresholdIntervals(group_mask, :)), [], 2) / 2;
nStdErrors = logErrors ./ logSigmas;

masked_groups = ternaryAvg.groups(group_mask);

%% Compare different kinds of errors

fig = figure;
fig.Units = 'inches';
fig.Position = [2 2 12 4];

plot_desc = {{'absolute errors', absErrors}, {'log errors', logErrors}, ...
    {'no. stdev.', nStdErrors}};

for i = 1:3
    subplot(1, 3, i);
    hold on;
    crt_plot_desc1 = plot_desc{i};
    crt_plot_desc2 = plot_desc{mod(i, 3) + 1};
    
    crt_data1 = crt_plot_desc1{2};
    crt_data2 = crt_plot_desc2{2};
    scatter(crt_data1, crt_data2, [], colorDict('blue'));
    plot(xlim, [0 0], 'k:', 'linewidth', 2);
    plot([0 0], ylim, 'k:', 'linewidth', 2);
    xlabel(crt_plot_desc1{1});
    ylabel(crt_plot_desc2{1});
    
    beautifygraph;
end

preparegraph;

% takeaway: abs and log errors are very strongly and monotonically related,
%           so we can use just one of them (I'll focus on log)
%           no. stdev. can be large when log error is small and vice versa,
%           so this needs to be looked into more carefully

%% Any differences between simple and mixed planes?

figure;

mask_mixed = cellfun(@(g) sum(g == ';') == 1, masked_groups);
mask_simple = ~mask_mixed;

hold on;
log_err_min = 1.1 * min(logErrors);
log_err_max = 1.1 * max(logErrors);
ci_n_std = 3;
fill(ci_n_std * [-1 -1 1 1 -1], ...
    [log_err_min log_err_max log_err_max log_err_min log_err_min], ...
    colorDict('light orange'), 'edgecolor', 'none');
scatter_opts = {'filled', 'markerfacealpha', 0.7, 'markeredgealpha', 0.7};
h_simple = scatter(nStdErrors(mask_simple), logErrors(mask_simple), [], ...
    colorDict('red'), scatter_opts{:}, 'displayname', 'simple planes');
h_mixed = scatter(nStdErrors(mask_mixed), logErrors(mask_mixed), [], ...
    colorDict('blue'), scatter_opts{:}, 'displayname', 'mixed planes');

plot(xlim, [0 0], 'k:', 'linewidth', 2);
plot([0 0], ylim, 'k:', 'linewidth', 2);

xlabel('no. stdev');
ylabel('log errors');

legend([h_simple h_mixed], 'location', 'southeast');
legend('boxoff');

ylim([log_err_min, log_err_max]);

beautifygraph;
preparegraph;

% takeaway: we mostly underestimate thresholds in simple planes
%           largest underestimates are actually within error bars