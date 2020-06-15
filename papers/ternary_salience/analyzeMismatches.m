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

% we know we won't get the first-order stats right, so get rid of those
% group_mask = ~strcmp(ternaryAvg.groups, 'A_1');

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

% takeaway: abs and log errors are very strongly and monotonically related,
%           so we can use just one of them (I'll focus on log)
%           no. stdev. can be large when log error is small and vice versa,
%           so this needs to be looked into more carefully

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

% takeaways: we mostly underestimate thresholds in simple planes
%            largest underestimates are actually within error bars

%% Look again at simple vs. mixed planes

mixed_directions = cell2mat(masked_directions(mask_mixed));
mixed_axis_dist = squeeze(std(reshape(mixed_directions', ...
    [3, 2, length(mixed_directions)]), [], 1));
on_axis_mixed0 = min(mixed_axis_dist, [], 1) < 1e-6;

on_axis_mixed = false(size(mask_mixed));
on_axis_mixed(mask_mixed) = on_axis_mixed0;

mask_mixed_not_on_axis = mask_mixed & ~on_axis_mixed;
mask_simple = ~mask_mixed;

% compare to previous plot after removing on-axis mixed-plane data
% figure;
% 
% hold on;
% log_err_min = 1.2 * min(logErrors);
% log_err_max = 1.2 * max(logErrors);
% ci_n_std = 3;
% fill(ci_n_std * [-1 -1 1 1 -1], ...
%     [log_err_min log_err_max log_err_max log_err_min log_err_min], ...
%     colorDict('light orange'), 'edgecolor', 'none');
% scatter_opts = {'filled', 'markerfacealpha', 0.7, 'markeredgealpha', 0.7};
% h_simple = scatter(nStdErrors(mask_simple), logErrors(mask_simple), [], ...
%     colorDict('red'), scatter_opts{:}, 'displayname', 'simple planes');
% h_mixed = scatter(nStdErrors(mask_mixed_not_on_axis), logErrors(mask_mixed_not_on_axis), [], ...
%     colorDict('blue'), scatter_opts{:}, 'displayname', 'mixed planes');
% 
% plot(xlim, [0 0], 'k:', 'linewidth', 2);
% plot([0 0], ylim, 'k:', 'linewidth', 2);
% 
% xlabel('no. stdev');
% ylabel('log errors');
% 
% legend([h_simple h_mixed], 'location', 'southeast');
% legend('boxoff');
% 
% ylim([log_err_min, log_err_max]);
% 
% beautifygraph;
% preparegraph;

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_simple = ksdensity(logErrors(mask_simple), kde_x);
kde_mixed = ksdensity(logErrors(mask_mixed_not_on_axis), kde_x);

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

disp(['simple vs. mixed planes: KS p = ' num2str(simple_mixed_ks, '%.2g') ...
    ', ranksum p = ' num2str(simple_mixed_ranksum, '%.2g')]);

% takeaways: simple-plane distribution is shifted towards negative values
%            simple-plane and mixed-plane distribution stat. sig. different
%               (according to both KS and ranksum tests)

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
kde_on_axis = ksdensity(logErrors(on_axis), kde_x);
kde_off_axis = ksdensity(logErrors(off_axis), kde_x);

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

disp(['on vs. off-axis: KS p = ' num2str(on_off_ks, '%.2g') ...
    ', ranksum p = ' num2str(on_off_ranksum, '%.2g')]);

% takeaways: on-axis and off-axis distributions look the same
%               (both KS and ranksum tests give p > 0.5)

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
kde_plus = ksdensity(logErrors(mask_plus), kde_x);
kde_minus = ksdensity(logErrors(mask_minus), kde_x);

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

disp(['plus vs. minus directions: KS p = ' num2str(plus_minus_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_ranksum, '%.2g')]);

% takeaways: there don't seem to be big differences between plus and
%            minus directions, though one could claim there are slightly
%            *larger* errors in the plus directions
%               (KS test yields p = 0.046; ranksum yields p = 0.54)

%% Plus (_11) vs. minus (_12) directions ignoring mixed plus-minus ones

group_indices = cellfun(@(s) s(find(s == '_') + 1), masked_groups, ...
    'UniformOutput', false);

mask_plus_exc = strcmp(group_indices, '11') | strcmp(group_indices, '1111');
mask_minus_exc = strcmp(group_indices, '12') | strcmp(group_indices, '1212');

% show KDE of signed errors
log_err_max_abs = max(abs(log_err_min), abs(log_err_max));
kde_x = linspace(-log_err_max_abs, log_err_max_abs, 250);
kde_plus_exc = ksdensity(logErrors(mask_plus_exc), kde_x);
kde_minus_exc = ksdensity(logErrors(mask_minus_exc), kde_x);

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

disp(['plus vs. pure minus directions: KS p = ' num2str(plus_minus_exc_ks, '%.2g') ...
    ', ranksum p = ' num2str(plus_minus_exc_ranksum, '%.2g')]);

% takeaways: there don't seem to be big differences between plus and
%            minus directions, though one could claim there are slightly
%            *larger* errors in the plus directions
%               (KS test yields p = 0.046; ranksum yields p = 0.54)

%% XXX separate out mixed planes like AB11/AC11 vs. AB11/AB12
