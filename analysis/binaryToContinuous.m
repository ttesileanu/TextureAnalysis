% Reanalyzing binary psychopphysics using continuous texture statistics.

%% Load psychophysics data

f = fopen(fullfile('data', 'btc_thresholds.csv'));
data0 = textscan(f, '%s', 'delimiter', '');
fclose(f);

% process the data a little

% keep only the stuff combined by hand for beta and theta
mask = true(1, 24);
%mask(end-5:end) = false;
mask(3:8) = false;
mask(11:18) = false;

% split at delimiters
group_type = strsplit(data0{1}{1}, ',', 'CollapseDelimiters', false);
group_type = group_type(2:end);
group_type = group_type(mask);

groups = strsplit(data0{1}{2}, ',', 'CollapseDelimiters', false);
groups = groups(2:end);
groups = groups(mask);

signs = strsplit(data0{1}{3}, ',', 'CollapseDelimiters', false);
signs = signs(2:end);
signs = signs(mask);

subjects = struct;
for i = 4:length(data0{1})
    crt_thresholds0 = strsplit(data0{1}{i}, ',', 'CollapseDelimiters', false);
    crt_thresholds = cellfun(@str2double, crt_thresholds0(2:end));
    subjects.(crt_thresholds0{1}) = crt_thresholds(mask);
end

% take geomean of positive and negative measurements
group_type_signmean = group_type(1:2:end);
groups_signmean = groups(1:2:end);
subjects_signmean = subjects;
fields = fieldnames(subjects_signmean);
for i = 1:length(fields)
    crt_thresholds = subjects_signmean.(fields{i});
    crt_thresholds_signmean = zeros(1, length(crt_thresholds)/2);
    for j = 1:length(crt_thresholds_signmean)
        j0 = 2*j - 1;
        crt_thresholds_signmean(j) = geomean(crt_thresholds(j0:j0+1));
    end
    subjects_signmean.(fields{i}) = crt_thresholds_signmean;
end

%% Load natural image stats

natural_stats0 = open(fullfile('save', 'natural_nosky_multiscale_60_AF_1_2_3_4_5_6.mat'));
blockAF_choice = 2;
natural_stats = natural_stats0.res{blockAF_choice};

% ann_stats = open('AnnNIstats.mat');
% natural_stats = struct('ev', ann_stats.NIstats.set(6).stats);

%% Predict thresholds

natural_cov = cov(natural_stats.ev);

% assuming threshold = variance^trafo_pow; or equivalently,
% sensitivity = variance^(-trafo_pow)
trafo_pow = -0.5;
predicted_thresholds = diag(natural_cov).^trafo_pow;
predicted_sensitivities = 1./predicted_thresholds;

%% Collect thresholds in groups

fields = fieldnames(subjects_signmean);
subjects_grouped_thresholds = struct;
for i = 1:length(fields)
    crt_thresholds = subjects_signmean.(fields{i});
    
    % gamma, beta-card, beta-diag, theta, alpha
    crt_grouped_thresholds = nan(1, 5);
    
    gamma_mask = strcmp(group_type_signmean, 'gamma');
    alpha_mask = strcmp(group_type_signmean, 'alpha');
    beta_card_mask = strcmp(group_type_signmean, 'beta-card') | strcmp(group_type_signmean, 'beta-card all');
    beta_diag_mask = strcmp(group_type_signmean, 'beta-diag') | strcmp(group_type_signmean, 'beta-diag all');
    theta_mask = strcmp(group_type_signmean, 'theta') | strcmp(group_type_signmean, 'theta all');
    
    if sum(gamma_mask) > 0
        crt_grouped_thresholds(1) = geomean(crt_thresholds(gamma_mask));
    end
    if sum(beta_card_mask) > 0
        crt_grouped_thresholds(2) = geomean(crt_thresholds(beta_card_mask));
    end
    if sum(beta_diag_mask) > 0
        crt_grouped_thresholds(3) = geomean(crt_thresholds(beta_diag_mask));
    end
    if sum(theta_mask) > 0
        crt_grouped_thresholds(4) = geomean(crt_thresholds(theta_mask));
    end
    if sum(alpha_mask) > 0
        crt_grouped_thresholds(5) = geomean(crt_thresholds(alpha_mask));
    end
    
    subjects_grouped_thresholds.(fields{i}) = ...
        struct('thresholds', crt_grouped_thresholds, 'sensitivities', 1./crt_grouped_thresholds);
end

%% Compare

figure;
hold on;

grouped_predicted_thresholds = nan(1, 5);
grouped_predicted_sensitivities = nan(1, 5);

% beta-card
grouped_predicted_thresholds(2) = geomean(predicted_thresholds(2:3));
grouped_predicted_sensitivities(2) = geomean(predicted_sensitivities(2:3));
% beta-diag
grouped_predicted_thresholds(3) = geomean(predicted_thresholds(4:5));
grouped_predicted_sensitivities(3) = geomean(predicted_sensitivities(4:5));
% theta
grouped_predicted_thresholds(4) = geomean(predicted_thresholds(6:9));
grouped_predicted_sensitivities(4) = geomean(predicted_sensitivities(6:9));
% alpha
grouped_predicted_thresholds(5) = predicted_thresholds(10);
grouped_predicted_sensitivities(5) = predicted_sensitivities(10);

fields = fieldnames(subjects_grouped_thresholds);
% colors = parula(length(fields));
colors = cool(length(fields));
for i = 1:length(fields)
    crt_sensitivities = subjects_grouped_thresholds.(fields{i}).sensitivities;
    crt_color = colors(i, :);
    smartscatter(grouped_predicted_sensitivities, crt_sensitivities, 250, crt_color, 'marker', '.');
    crt_mask = isfinite(crt_sensitivities(:)) & isfinite(grouped_predicted_sensitivities(:));
    drawfitline(grouped_predicted_sensitivities(crt_mask), crt_sensitivities(crt_mask), 'style', ...
        {'linewidth', 1, 'color', crt_color}, 'legend', false, ...
        'showci', false);
%    plot(predicted_thresholds, crt_thresholds, '.', 'color', crt_color, 'markersize', 25);
    disp(['Correlation ' fields{i} ': ' ...
        num2str(corr(flatten(crt_sensitivities(crt_mask)), flatten(grouped_predicted_sensitivities(crt_mask))), '%.2f')]);
end
xlabel('Predicted sensitivities');
ylabel('Measured sensitivities');

beautifygraph;
preparegraph;

figure;

c_labels = {'\gamma', '\beta_|', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{|-}', '\theta_{-|}', '\theta_{\_|}', '\theta_{|\_}', '\alpha'};

% find optimal scaling
all_measured_thresholds = zeros(length(fields), length(grouped_predicted_thresholds));
all_measured_sensitivities = zeros(length(fields), length(grouped_predicted_thresholds));
for i = 1:length(fields)
    all_measured_thresholds(i, :) = subjects_grouped_thresholds.(fields{i}).thresholds;
    all_measured_sensitivities(i, :) = subjects_grouped_thresholds.(fields{i}).sensitivities;
end

% need to ignore gamma
average_measured_thresholds = geomean(all_measured_thresholds, 1);
average_measured_sensitivities = geomean(all_measured_sensitivities, 1);
scaling_thresholds = ...
    dot(grouped_predicted_thresholds(2:end), average_measured_thresholds(2:end)) / ...
    norm(grouped_predicted_thresholds(2:end))^2;
scaling_sensitivities = ...
    dot(grouped_predicted_sensitivities(2:end), average_measured_sensitivities(2:end)) / ...
     norm(grouped_predicted_sensitivities(2:end))^2;

hold on;
% locations of the 9 statistics (matching c_labels except gamma) on the x axis
all_x_positions = [1:2 5:6 9:12 15];
plot(all_x_positions, predicted_sensitivities(2:end)*scaling_sensitivities, 'bo', ...
    'markersize', 15);

measurement_x_positions = [3, 7, 13, 16];
for i = 1:length(fields)
    plot(measurement_x_positions, subjects_grouped_thresholds.(fields{i}).sensitivities(2:end), ...
        's', 'color', colors(i, :), 'markersize', 15);
end
ylabel('Sensitivity');

xlim([0 17]);
ylim([0.5 4.5]);

beautifygraph;

set(gca, 'xtick', all_x_positions, 'xticklabel', c_labels(2:end), 'xminortick', 'off');

preparegraph;