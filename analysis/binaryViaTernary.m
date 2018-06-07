% Analyzing binary PP data using ternary stats.

%% Load binary psychophysics data

[binary, binary_raw] = loadBinaryPP(fullfile('data', 'PPandNIstats.mat'), 'ninterp', 3);
binary_by_subject = loadBinaryPP(fullfile('data', 'PPandNIstats.mat'), ...
    'subjects', '*', 'ninterp', 3);

%% Generate binary to ternary mapping

% This should be doable analytically. Let's see.
% 
% \gamma:
%   A_1:
%       p(A_1 == 0) = n_-/(n_+ + n_-) = 1/2*(1 - \gamma)
%       p(A_1 == 1) = 0
%       p(A_1 == 2) = n_+/(n_+ + n_-) = 1/2*(1 + \gamma)
%   A_2:
%       p(A_2 == 0) = P(A == 0) = 1/2*(1 - \gamma)
%       p(A_2 == 1) = P(A == 2) = 1/2*(1 + \gamma)
%       p(A_2 == 2) = P(A == 1) = 0
%   p(AB_11 == 0) = p(A == 0)*p(A == 0) + p(A == 1)*p(A == 2) + p(A == 2)*p(A == 1)
%                 = 1/4*(1 - \gamma)^2
%   p(AB_11 == 1) = p(A == 0)*p(A == 1) + p(A == 1)*p(A == 0) + p(A == 2)*p(A == 2)
%                 = 1/4*(1 + \gamma)^2
%   p(AB_11 == 2) = p(A == 0)*p(A == 2) + p(A == 2)*p(A == 0) + p(A == 1)*p(A == 1)
%                 = 1/2*(1 - \gamma^2)
%
% \beta_--:
%   p(AB_11 == 0) = 

nlocs = 8;
nsamples_binary = 6;
patch_size = 64;

t0 = tic;
set(groot, 'defaultFigureVisible', 'off');
[binary_mapping, binary_map_details] = generateTextureMapping(....
    @(patch) expandev(flatten(getnth(2, @processBlock, patch, 3))', 3), ...
    binary.groups, binary.directions, ...
    'nlocs', nlocs, 'nsamples', nsamples_binary, 'patchsize', patch_size);
set(groot, 'defaultFigureVisible', 'on');
disp(['Binary-to-ternary map generation took ' num2str(toc(t0), '%.1f') ' seconds.']);

%% Save binary to ternary mapping

save(fullfile('save', 'binary_to_ternary_map.mat'), 'binary_mapping', 'binary_map_details', ...
    'nlocs', 'nsamples_binary', 'patch_size', 'binary');

%% Load binary to ternary mapping

load(fullfile('save', 'binary_to_ternary_map.mat'));

%% Figure out gains from ternary statistics

% natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_contrastadapt_with_focus.mat'));
natural_stats0 = open(fullfile('save', 'natural_nosky_ternary_with_focus.mat'));
natural_stats = natural_stats0.res{1};

restrict_to_focus = true;

% restrict to in-focus patches, if possible
if restrict_to_focus && isfield(natural_stats, 'focus')
    disp('Restricting to in-focus patches.');
    mask = (natural_stats.focus.clusterIds == natural_stats.focus.focusCluster);
    fields = {'objIds', 'ev', 'pxPerPatch', 'patchLocations', 'patchLocationsOrig', 'imgIds'};
    for i = 1:length(fields)
        natural_stats.(fields{i}) = natural_stats.(fields{i})(mask, :);
    end
    natural_stats.covM = cov(natural_stats.ev);
end

% expand natural stats to contain all 99 probabilities
ni_ev_full = expandev(natural_stats.ev, 3);

% load ternary psychophysics, used only to fix scale for gain matrix
ternary_avg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% calculate unscaled gain matrix
[gain0, ~, pred_details] = getTernaryPredictions(ni_ev_full, ...
    ternary_avg, 2.0*eye(size(ni_ev_full, 2)), eye(size(ni_ev_full, 2)), 3e-5, ...
    'fitscale_opts', {'exclude', ...
    strcmp(ternary_avg.groups, 'A_1') & cellfun(@length, ternary_avg.groups) > 6}, ...
    'effcode_opts', {'varsal', true});

% scale the gain matrix appropriately
gain = gain0 / pred_details.acoeff;

%% Check that the mapping makes sense

planes = sortgroups(unique(ternary_avg.groups(~ternary_avg.multi)));
to_display = 'BC_1_1';
crt_idx = find(strcmp(binary.groups, to_display), 1);
crt_mapping = binary_mapping{crt_idx}.function;

locs_to_display = linspace(0, 1, 8);

plotter = MatrixPlotter(length(planes));
while plotter.next
    i = plotter.index;
    
    hold on;
    drawTernaryTriangle;
    
    crt_group = planes{i};
    crt_dirs = ternaryproject(crt_mapping(locs_to_display)', crt_group)';
    crt_projections = ternary3to2(crt_dirs);
%     scatter(crt_projections(:, 1), crt_projections(:, 2), [], locs_to_display, 'filled');
    smartscatter(crt_projections(:, 1), crt_projections(:, 2), ...
        'color', locs_to_display, 'density', false);
    
    axis equal;
    
    title(crt_group);
    
    set(gca, 'box', 'on', 'xminortick', 'on', 'yminortick', 'on', 'linewidth', 1);
end

set(gcf, 'Name', to_display);

%% Check that the mapping makes sense (mixed groups)

planes = sortgroups(unique(ternary_avg.groups(ternary_avg.multi)));
% planes = {'AB_1_1[0];AC_1_1[0]', 'AB_1_1[1];AC_1_1[1]', 'AB_1_1[2];AC_1_1[2]', ...
%     'AB_1_1[0];AC_1_1[1]', 'AB_1_1[1];AC_1_1[2]', 'AB_1_1[2];AC_1_1[0]', ...
%     'AB_1_1[0];AC_1_1[2]', 'AB_1_1[1];AC_1_1[0]', 'AB_1_1[2];AC_1_1[1]'};
to_display = 'AB_1_1;AC_1_1';
crt_idxs = find(strcmp(binary.groups, to_display));

locs_to_display = linspace(0, 1, 8);

plotter = MatrixPlotter(length(planes));
while plotter.next
    i = plotter.index;
    
    crt_group = planes{i};
    crt_groups = strtrim(strsplit(crt_group, ';'));
    
    hold on;
    drawTernaryMixedBackground(crt_groups{:});
    
    for k = 1:length(crt_idxs)
        crt_idx = crt_idxs(k);
        crt_mapping = binary_mapping{crt_idx}.function;
        
        crt_dirs = ternaryproject(crt_mapping(locs_to_display)', crt_group)';
        crt_projections = (3/2)*crt_dirs - 1/2;
        smartscatter(crt_projections(:, 1), crt_projections(:, 2), ...
            'color', locs_to_display, 'density', false);
    end
    
    axis equal;
    
    title(crt_group);
    
    set(gca, 'box', 'on', 'xminortick', 'on', 'yminortick', 'on', 'linewidth', 1);
end

set(gcf, 'Name', to_display);

%% Predict binary thresholds

% XXX for now taking numerical derivative
% XXX I think all the dependencies are linear though, and could be
% calculated analytically
num_der_dist = 0.01;
binary_predictions0 = gainsToThresholds(gain, cellfun(@(s) ...
    (s.function(num_der_dist) - s.function(0)) / num_der_dist, binary_mapping, ...
    'uniform', false));

% rescale binary predictions to match data as much as possible
[binary_acoeff, binary_predictions] = fitscale(binary_predictions0, binary.thresholds);
% binary_predictions = binary_predictions0;

binary_predictions_by_subject = matchThresholds(binary_predictions, ...
    binary.groups, binary.directions, ...
    binary_by_subject.groups, binary_by_subject.directions);

%% Check single-group binary predictions

figure;
plotBinaryThresholds(binary_predictions_by_subject, binary_by_subject);

%% Check multi-group binary predictions

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 12 12];
hold on;

plotBinaryEllipses(binary_predictions_by_subject, binary_by_subject, ...
    'predellipses', false);

preparegraph;