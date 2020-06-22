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

% keep only second-order groups (simple or mixed)
group_mask = cellfun(@(s) length(s) == 6 || sum(s == ';') == 1, ternaryAvg.groups);

%% Find absolute and log errors, and log errors divided by uncertainties

absErrors = predictions.thresholds(group_mask) - ternaryAvg.thresholds(group_mask);
logErrors = log(predictions.thresholds(group_mask)) - ...
    log(ternaryAvg.thresholds(group_mask));
logSigmas = diff(log(ternaryAvg.thresholdIntervals(group_mask, :)), [], 2) / 2;
nStdErrors = logErrors ./ logSigmas;

meanThresholds = 0.5 * (predictions.thresholds(group_mask) + ...
    ternaryAvg.thresholds(group_mask));
errorRatios = absErrors ./ meanThresholds;

masked_groups = ternaryAvg.groups(group_mask);
masked_directions = ternaryAvg.directions(group_mask);
mask_mixed = cellfun(@(g) sum(g == ';') == 1, masked_groups);

%% Compare different kinds of errors

fig = figure;
fig.Units = 'inches';
fig.Position = [2 2 16 4];

plot_desc = {{'absolute errors', absErrors}, {'log errors', logErrors}, ...
    {'no. stdev.', nStdErrors}};

for i = 1:3
    subplot(1, 4, i);
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

subplot(1, 4, 4);
hold on;
scatter(logErrors, errorRatios, [], colorDict('blue'));
plot(xlim, xlim, 'k--', 'linewidth' ,1);
plot(xlim, [0 0], 'k:', 'linewidth', 2);
plot([0 0], ylim, 'k:', 'linewidth', 2);
xlabel('log errors');
ylabel('symm. frac. errors');
beautifygraph;

preparegraph;

% takeaway: * abs and log errors are very strongly and monotonically
%             related, so we can use just one of them (I'll focus on log)
%           * no. stdev. can be large when log error is small and vice
%             versa, so this might need to be looked into more carefully
%           * symmetric percentage errors are essentially identical to log
%             errors (log errors are always larger, but maximum ratio
%             between the two is just 1.034)

%% Any differences between simple and mixed planes?

figure;

mask_simple = ~mask_mixed;

hold on;
log_err_min = 1.2 * min(logErrors);
log_err_max = 1.2 * max(logErrors);
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

% takeaways: * we mostly underestimate thresholds in simple planes
%            * largest underestimates are mostly within error bars

%% Store all the results in a map

all_results = containers.Map;

%% Look again at simple vs. mixed planes

mixed_directions = cell2mat(masked_directions(mask_mixed));
mixed_axis_dist = squeeze(std(reshape(mixed_directions', ...
    [3, 2, length(mixed_directions)]), [], 1));
on_axis_mixed0 = min(mixed_axis_dist, [], 1) < 1e-6;

on_axis_mixed = false(size(mask_mixed));
on_axis_mixed(mask_mixed) = on_axis_mixed0;

mask_mixed_not_on_axis = mask_mixed & ~on_axis_mixed;
mask_simple = ~mask_mixed;

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.07};
kde_simple = ksdensity(logErrors(mask_simple), kde_x, kde_opts{:});
kde_mixed = ksdensity(logErrors(mask_mixed_not_on_axis), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_simple(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_mixed(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_simple, 'color', colorDict('red'));
plot(kde_x, kde_mixed, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'simple planes', 'mixed planes'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, simple_mixed_ks] = kstest2(logErrors(mask_simple), logErrors(mask_mixed_not_on_axis));
simple_mixed_ranksum = ranksum(logErrors(mask_simple), logErrors(mask_mixed_not_on_axis));
[~, simple_mixed_abs_ks] = kstest2(abs(logErrors(mask_simple)), ...
    abs(logErrors(mask_mixed_not_on_axis)));
simple_mixed_abs_ranksum = ranksum(abs(logErrors(mask_simple)), ...
    abs(logErrors(mask_mixed_not_on_axis)));

simple_mixed_res = struct;
simple_mixed_res.kde1 = kde_simple;
simple_mixed_res.kde2 = kde_mixed;
simple_mixed_res.label1 = 'simple';
simple_mixed_res.label2 = 'mixed';
simple_mixed_res.errors1 = logErrors(mask_simple);
simple_mixed_res.errors2 = logErrors(mask_mixed);
simple_mixed_res.p_ks = simple_mixed_ks;
simple_mixed_res.p_ranksum = simple_mixed_ranksum;
simple_mixed_res.p_ks_abs = simple_mixed_abs_ks;
simple_mixed_res.p_ranksum_abs = simple_mixed_abs_ranksum;

all_results('simple_vs_mixed') = simple_mixed_res;

disp(['simple vs. mixed planes: KS p = ' num2str(simple_mixed_ks, '%.2g') ...
    ', ranksum p = ' num2str(simple_mixed_ranksum, '%.2g')]);

% takeaways: * simple-plane distribution is shifted towards negative values
%            * simple-plane and mixed-plane distribution stat. sig.
%              different (according to both KS and ranksum tests)

%% On-axis vs. off-axis

simple_directions = cell2mat(masked_directions(mask_simple));
on_axis_simple0 = abs(max(abs(simple_directions), [], 2) - 1) < 1e-6;

on_axis_simple = false(size(mask_simple));
on_axis_simple(mask_simple) = on_axis_simple0;

on_axis = false(size(mask_simple));
on_axis(on_axis_simple) = true;
on_axis(on_axis_mixed) = true;

off_axis = ~on_axis;

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.07};
kde_on_axis = ksdensity(logErrors(on_axis), kde_x, kde_opts{:});
kde_off_axis = ksdensity(logErrors(off_axis), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_on_axis(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_off_axis(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_on_axis, 'color', colorDict('red'));
plot(kde_x, kde_off_axis, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'on axis', 'off axis'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, on_off_ks] = kstest2(logErrors(on_axis), logErrors(off_axis));
on_off_ranksum = ranksum(logErrors(on_axis), logErrors(off_axis));
[~, on_off_abs_ks] = kstest2(abs(logErrors(on_axis)), abs(logErrors(off_axis)));
on_off_abs_ranksum = ranksum(abs(logErrors(on_axis)), abs(logErrors(off_axis)));

on_off_res = struct;
on_off_res.kde1 = kde_on_axis;
on_off_res.kde2 = kde_off_axis;
on_off_res.label1 = 'on-axis';
on_off_res.label2 = 'off-axis';
on_off_res.errors1 = logErrors(on_axis);
on_off_res.errors2 = logErrors(off_axis);
on_off_res.p_ks = on_off_ks;
on_off_res.p_ranksum = on_off_ranksum;
on_off_res.p_ks_abs = on_off_abs_ks;
on_off_res.p_ranksum_abs = on_off_abs_ranksum;

all_results('on_vs_off') = on_off_res;

disp(['on vs. off-axis: KS p = ' num2str(on_off_ks, '%.2g') ...
    ', ranksum p = ' num2str(on_off_ranksum, '%.2g')]);

% takeaways: * on-axis and off-axis distributions look the same
%                (both KS and ranksum tests give p > 0.5)

%% Plus (_11) vs. minus (_12) directions

% for simple planes, simply find all the digits, and see whether we have
% one unique digit (which will always be 1) or two (which will be 1 and 2)
% the former = plus direction, latter = minus direction
mask_plus = false(size(mask_simple));
mask_plus(mask_simple) = cellfun(@(s) ...
    length(unique(s(isstrprop(s, 'digit')))), ...
    masked_groups(mask_simple)) == 1;

% for mixed planes this is complicated by the numbers is square brackets
% to get around that, find all the underscores, and look at characters one
% past them --> these are the digits identifying the directions
% now we can check whether there's only one unique such digit ( --> plus)
% or two ( --> minus direction)
mask_plus(mask_mixed) = cellfun(@(s) ...
    length(unique(s(find(s == '_') + 1))), ...
    masked_groups(mask_mixed)) == 1;

mask_minus = ~mask_plus;

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.07};
kde_plus = ksdensity(logErrors(mask_plus), kde_x, kde_opts{:});
kde_minus = ksdensity(logErrors(mask_minus), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_plus(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_minus(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_plus, 'color', colorDict('red'));
plot(kde_x, kde_minus, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'''plus'' (11) directions', '''minus'' (12) directions'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, plus_minus_ks] = kstest2(logErrors(mask_plus), logErrors(mask_minus));
plus_minus_ranksum = ranksum(logErrors(mask_plus), logErrors(mask_minus));
[~, plus_minus_abs_ks] = kstest2(abs(logErrors(mask_plus)), ...
    abs(logErrors(mask_minus)));
plus_minus_abs_ranksum = ranksum(abs(logErrors(mask_plus)), ...
    abs(logErrors(mask_minus)));

plus_minus_res = struct;
plus_minus_res.kde1 = kde_plus;
plus_minus_res.kde2 = kde_minus;
plus_minus_res.label1 = '''sum''';
plus_minus_res.label2 = '''difference''';
plus_minus_res.errors1 = logErrors(mask_plus);
plus_minus_res.errors2 = logErrors(mask_minus);
plus_minus_res.p_ks = plus_minus_ks;
plus_minus_res.p_ranksum = plus_minus_ranksum;
plus_minus_res.p_ks_abs = plus_minus_abs_ks;
plus_minus_res.p_ranksum_abs = plus_minus_abs_ranksum;

all_results('plus_vs_minus') = plus_minus_res;

disp(['plus vs. minus directions: KS p = ' num2str(plus_minus_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_ranksum, '%.2g')]);

% takeaways: * there don't seem to be big differences between plus and
%              minus directions, though one could claim there are slightly
%              *larger* errors in the plus directions
%                (KS test yields p = 0.046; ranksum yields p = 0.54)

%% Plus (_11) vs. minus (_12) directions ignoring mixed plus-minus ones

group_indices = cellfun(@(s) s(find(s == '_') + 1), masked_groups, ...
    'UniformOutput', false);

mask_plus_exc = strcmp(group_indices, '11') | strcmp(group_indices, '1111');
mask_minus_exc = strcmp(group_indices, '12') | strcmp(group_indices, '1212');

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.07};
kde_plus_exc = ksdensity(logErrors(mask_plus_exc), kde_x, kde_opts{:});
kde_minus_exc = ksdensity(logErrors(mask_minus_exc), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_plus_exc(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_minus_exc(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_plus_exc, 'color', colorDict('red'));
plot(kde_x, kde_minus_exc, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'''plus'' (11) directions', '''pure'' ''minus'' (12) directions'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, plus_minus_exc_ks] = kstest2(logErrors(mask_plus_exc), logErrors(mask_minus_exc));
plus_minus_exc_ranksum = ranksum(logErrors(mask_plus_exc), logErrors(mask_minus_exc));
[~, plus_minus_exc_abs_ks] = kstest2(abs(logErrors(mask_plus_exc)), abs(logErrors(mask_minus_exc)));
plus_minus_exc_abs_ranksum = ranksum(abs(logErrors(mask_plus_exc)), abs(logErrors(mask_minus_exc)));

plus_minus_exc_res = struct;
plus_minus_exc_res.kde1 = kde_plus_exc;
plus_minus_exc_res.kde2 = kde_minus_exc;
plus_minus_exc_res.label1 = '''sum''';
plus_minus_exc_res.label2 = '''difference''';
plus_minus_exc_res.errors1 = logErrors(mask_plus_exc);
plus_minus_exc_res.errors2 = logErrors(mask_minus_exc);
plus_minus_exc_res.p_ks = plus_minus_exc_ks;
plus_minus_exc_res.p_ranksum = plus_minus_exc_ranksum;
plus_minus_exc_res.p_ks_abs = plus_minus_exc_abs_ks;
plus_minus_exc_res.p_ranksum_abs = plus_minus_exc_abs_ranksum;

all_results('plus_vs_minus_exc') = plus_minus_exc_res;

disp(['plus vs. pure minus directions: KS p = ' num2str(plus_minus_exc_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_exc_ranksum, '%.2g')]);

% takeaways: * there don't seem to be big differences between plus and
%              minus directions, though one could claim there are slightly
%              *larger* errors in the plus directions
%                (KS test yields p = 0.046; ranksum yields p = 0.54)

%% Plus (_11) vs. minus (_12) directions in simple planes

group_indices = cellfun(@(s) s(find(s == '_') + 1), masked_groups, ...
    'UniformOutput', false);

mask_plus_simple = strcmp(group_indices, '11');
mask_minus_simple = strcmp(group_indices, '12');

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.12};
kde_plus_simple = ksdensity(logErrors(mask_plus_simple), kde_x, kde_opts{:});
kde_minus_simple = ksdensity(logErrors(mask_minus_simple), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_plus_simple(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_minus_simple(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_plus_simple, 'color', colorDict('red'));
plot(kde_x, kde_minus_simple, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'simple ''plus'' (11) directions', 'simple ''minus'' (12) directions'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, plus_minus_simple_ks] = kstest2(logErrors(mask_plus_simple), logErrors(mask_minus_simple));
plus_minus_simple_ranksum = ranksum(logErrors(mask_plus_simple), logErrors(mask_minus_simple));
[~, plus_minus_simple_abs_ks] = kstest2(abs(logErrors(mask_plus_simple)), ...
    abs(logErrors(mask_minus_simple)));
plus_minus_simple_abs_ranksum = ranksum(abs(logErrors(mask_plus_simple)), ...
    abs(logErrors(mask_minus_simple)));

plus_minus_simple_res = struct;
plus_minus_simple_res.kde1 = kde_plus_simple;
plus_minus_simple_res.kde2 = kde_minus_simple;
plus_minus_simple_res.label1 = 'simple ''sum''';
plus_minus_simple_res.label2 = 'simple ''difference''';
plus_minus_simple_res.errors1 = logErrors(mask_plus_simple);
plus_minus_simple_res.errors2 = logErrors(mask_minus_simple);
plus_minus_simple_res.p_ks = plus_minus_simple_ks;
plus_minus_simple_res.p_ranksum = plus_minus_simple_ranksum;
plus_minus_simple_res.p_ks_abs = plus_minus_simple_abs_ks;
plus_minus_simple_res.p_ranksum_abs = plus_minus_simple_abs_ranksum;

all_results('plus_vs_minus_simple') = plus_minus_simple_res;

disp(['simple-plane plus vs. minus directions: KS p = ' num2str(plus_minus_simple_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_simple_ranksum, '%.2g')]);

% takeaways: * there don't seem to be big differences between simple-plane
%              plus and minus directions, either
%                (KS test yields p = 0.22; ranksum yields p = 0.72)

%% Plus (_11) vs. minus (_12) directions in mixed planes

group_indices = cellfun(@(s) s(find(s == '_') + 1), masked_groups, ...
    'UniformOutput', false);

mask_plus_mixed = strcmp(group_indices, '1111');
mask_minus_mixed = strcmp(group_indices, '1212');

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_opts = {'Bandwidth', 0.07};
kde_plus_mixed = ksdensity(logErrors(mask_plus_mixed), kde_x, kde_opts{:});
kde_minus_mixed = ksdensity(logErrors(mask_minus_mixed), kde_x, kde_opts{:});

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_plus_mixed(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_minus_mixed(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_plus_mixed, 'color', colorDict('red'));
plot(kde_x, kde_minus_mixed, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'mixed ''plus'' (11) directions', 'mixed ''pure'' ''minus'' (12) directions'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, plus_minus_mixed_ks] = kstest2(logErrors(mask_plus_mixed), logErrors(mask_minus_mixed));
plus_minus_mixed_ranksum = ranksum(logErrors(mask_plus_mixed), logErrors(mask_minus_mixed));
[~, plus_minus_mixed_abs_ks] = kstest2(abs(logErrors(mask_plus_mixed)), abs(logErrors(mask_minus_mixed)));
plus_minus_mixed_abs_ranksum = ranksum(abs(logErrors(mask_plus_mixed)), abs(logErrors(mask_minus_mixed)));

plus_minus_mixed_res = struct;
plus_minus_mixed_res.kde1 = kde_plus_mixed;
plus_minus_mixed_res.kde2 = kde_minus_mixed;
plus_minus_mixed_res.label1 = 'mixed ''sum''';
plus_minus_mixed_res.label2 = 'mixed ''difference''';
plus_minus_mixed_res.errors1 = logErrors(mask_plus_mixed);
plus_minus_mixed_res.errors2 = logErrors(mask_minus_mixed);
plus_minus_mixed_res.p_ks = plus_minus_mixed_ks;
plus_minus_mixed_res.p_ranksum = plus_minus_mixed_ranksum;
plus_minus_mixed_res.p_ks_abs = plus_minus_mixed_abs_ks;
plus_minus_mixed_res.p_ranksum_abs = plus_minus_mixed_abs_ranksum;

all_results('plus_vs_minus_mixed') = plus_minus_mixed_res;

disp(['mixed-plane plus vs. pure minus directions: KS p = ' num2str(plus_minus_mixed_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_mixed_ranksum, '%.2g')]);

% takeaways: * errors seem to be larger in mixed plus directions
%                (KS p = 0.00032, ranksum p = 0.23)

%% 2d mixed planes (like AB/AC) vs. 1d mixed planes (like AB/AB)

mask_mixed_2d = cellfun(...
    @(s) sum(s == ';') == 1 && length(unique(s(isletter(s)))) == 3, ...
    masked_groups ...
);
mask_mixed_1d = cellfun(...
    @(s) sum(s == ';') == 1 && length(unique(s(isletter(s)))) == 2, ...
    masked_groups ...
);

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_mixed_2d = ksdensity(logErrors(mask_mixed_2d), kde_x);
kde_mixed_1d = ksdensity(logErrors(mask_mixed_1d), kde_x);

figure;
hold on;

fill([kde_x(:) ; flipud(kde_x(:))], [kde_mixed_2d(:) ; zeros(size(kde_x(:)))], ...
    colorDict('red'), 'facealpha', 0.3);
fill([kde_x(:) ; flipud(kde_x(:))], [kde_mixed_1d(:) ; zeros(size(kde_x(:)))], ...
    colorDict('blue'), 'facealpha', 0.3);

plot(kde_x, kde_mixed_2d, 'color', colorDict('red'));
plot(kde_x, kde_mixed_1d, 'color', colorDict('blue'));

xlabel('log errors');
ylabel('pdf');

legend({'2d mixed planes', '1d mixed planes'});
legend('boxoff');

beautifygraph;
preparegraph;

[~, mixed_2d_1d_ks] = kstest2(logErrors(mask_mixed_2d), logErrors(mask_mixed_1d));
mixed_2d_1d_ranksum = ranksum(logErrors(mask_mixed_2d), logErrors(mask_mixed_1d));
[~, mixed_2d_1d_abs_ks] = kstest2(abs(logErrors(mask_mixed_2d)), ...
    abs(logErrors(mask_mixed_1d)));
mixed_2d_1d_abs_ranksum = ranksum(abs(logErrors(mask_mixed_2d)), ...
    abs(logErrors(mask_mixed_1d)));

mixed_2d_1d_res = struct;
mixed_2d_1d_res.kde1 = kde_mixed_1d;
mixed_2d_1d_res.kde2 = kde_mixed_2d;
mixed_2d_1d_res.label1 = '''1d'' mixed';
mixed_2d_1d_res.label2 = '''2d'' mixed';
mixed_2d_1d_res.errors1 = logErrors(mask_mixed_1d);
mixed_2d_1d_res.errors2 = logErrors(mask_mixed_2d);
mixed_2d_1d_res.p_ks = mixed_2d_1d_ks;
mixed_2d_1d_res.p_ranksum = mixed_2d_1d_ranksum;
mixed_2d_1d_res.p_ks_abs = mixed_2d_1d_abs_ks;
mixed_2d_1d_res.p_ranksum_abs = mixed_2d_1d_abs_ranksum;

all_results('mixed_1d_vs_2d') = mixed_2d_1d_res;

disp(['2d vs. 1d mixed planes: KS p = ' num2str(mixed_2d_1d_ks, '%.2g') ...
    ', ranksum p = ' num2str(mixed_2d_1d_ranksum, '%.2g')]);

% takeaways: * no big differences between 1d and 2d mixed-plane errors
%                (KS test yields p = 0.25, ranksum yields p = 0.96)
%            * one might argue that large negative prediction errors in the
%              1d planes are more common than in the 2d planes (by eye)
%            * absolute log errors are somewhat larger in 1d mixed planes
%                (KS and ranksum tests both yield p = 0.01)

%% Make figure

plot_choices = {...
    'simple_vs_mixed', 'plus_vs_minus_exc', 'on_vs_off', 'mixed_1d_vs_2d' ...
};

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.2 2.3];

n_plots = length(plot_choices);

for i = 1:n_plots
    crt_name = plot_choices{i};
    crt_res = all_results(crt_name);
    
    crt_kde1 = crt_res.kde1;
    crt_kde2 = crt_res.kde2;
    
    subplot(2, 2, i);
    hold on;
    
    fill([kde_x(:) ; flipud(kde_x(:))], [crt_kde1(:) ; zeros(size(kde_x(:)))], ...
        colorDict('red'), 'facealpha', 0.3);
    fill([kde_x(:) ; flipud(kde_x(:))], [crt_kde2(:) ; zeros(size(kde_x(:)))], ...
        colorDict('blue'), 'facealpha', 0.3);

    plot(kde_x, crt_kde1, 'color', colorDict('red'));
    plot(kde_x, crt_kde2, 'color', colorDict('blue'));
    
%     xlabel('log errors');
    ylabel('pdf');

%     legend({crt_res.label1, crt_res.label2});
%     legend('boxoff');

    crt_ax = gca;
    crt_ax.Units = 'points';
    crt_ax_pos_pt = crt_ax.Position;
    
    crt_xlim = xlim;
    crt_ylim = ylim;
    
    crt_means = [nanmean(crt_res.errors1) nanmean(crt_res.errors2)];
    crt_medians = [nanmedian(crt_res.errors1) nanmedian(crt_res.errors2)];
    
    crt_color_order = {'red', 'blue'};
    for k = 1:2
        crt_color = colorDict(crt_color_order{k});
        crt_mean = crt_means(k);
        crt_median = crt_medians(k);
        
%         plot([crt_mean crt_mean], crt_ylim, 'color', crt_color);
        plot([crt_median crt_median], crt_ylim, 'color', crt_color);
    end
    
    stem(crt_res.errors1, crt_ylim(2) * 0.1 * ones(size(crt_res.errors1)), ...
        '-','Marker','none', 'color', colorDict('red'));
    stem(crt_res.errors2, -crt_ylim(2) * 0.1 * ones(size(crt_res.errors2)), ...
        '-','Marker','none', 'color', colorDict('blue'));
    
    crt_ax.Clipping = 'off';
    
%     ylim([-0.1 * crt_ylim(2), crt_ylim(2)]);
    ylim([0, crt_ylim(2)]);
    
    pt_to_data_y = diff(crt_ylim) / crt_ax_pos_pt(4);
    text(...
        crt_xlim(2), crt_ylim(2), crt_res.label1, ...
        'color', colorDict('red'), 'fontsize', 8, 'fontname', 'helvetica', ...
        'horizontalalignment', 'right', 'verticalalignment', 'top' ...
    );
    text(...
        crt_xlim(2), crt_ylim(2) - 10 * pt_to_data_y, crt_res.label2, ...
        'color', colorDict('blue'), 'fontsize', 8, 'fontname', 'helvetica', ...
        'horizontalalignment', 'right', 'verticalalignment', 'top' ...
    );

    beautifygraph;
    
    set(crt_ax, 'color','none');
end

preparegraph;

set(fig, 'color', 'none');
safePrint(fullfile('figs', 'draft', 'mismatchPlots.pdf'));

%% Make summary table

res_header = {'comparison', 'p_KS', 'p_ranksum', 'p_KS_abs', ...
    'p_ranksum_abs', 'median1', 'median2', 'sem1', 'sem2', ...
    'median_abs1', 'median_abs2', 'sem_abs1', 'sem_abs2'};
row_names = {};
res_array = [];
sem = @(v) nanstd(v) / sqrt(sum(isfinite(v)));
for i = 1:n_plots
    crt_name = plot_choices{i};
    crt_res = all_results(crt_name);
    
    row_names = [row_names crt_name]; %#ok<AGROW>
    crt_row = [crt_res.p_ks, crt_res.p_ranksum, ...
        crt_res.p_ks_abs, crt_res.p_ranksum_abs, ...
        nanmedian(crt_res.errors1), nanmedian(crt_res.errors2), ...
        sem(crt_res.errors1), sem(crt_res.errors2), ...
        nanmedian(abs(crt_res.errors1)), nanmedian(abs(crt_res.errors2)), ...
        sem(abs(crt_res.errors1)), sem(abs(crt_res.errors2)), ...
        ];
    res_array = [res_array ; crt_row]; %#ok<AGROW>
end

res_table = table(row_names', res_array(:, 1), res_array(:, 2), ...
    res_array(:, 3), res_array(:, 4), res_array(:, 5), res_array(:, 6), ...
    res_array(:, 7), res_array(:, 8), res_array(:, 9), res_array(:, 10), ...
    res_array(:, 11), res_array(:, 12), 'variablenames', res_header);
res_table %#ok<NOPTS>
