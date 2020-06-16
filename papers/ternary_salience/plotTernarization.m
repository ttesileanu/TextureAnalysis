% plot results of changing ternarization procedure

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   symmetrizePP
%       Set to `true` to take the average between each psychophysics
%       measurement and the measurement in the opposite texture direction.
%       This effectively forces the measurements to be centered at the
%       origin.
%   plotType
%       This can be 'violin' or 'jitter'.
%   diffType
%       Type of error to calculate. This can be
%        'spe'
%           Symmetric percentage error -- this is the ratio between the
%           difference and the mean of the measured and predicted values.
%        'log'
%           Difference between log measurement and log prediction.

setdefault('symmetrizePP', false);
setdefault('plotType', 'violin');
setdefault('diffType', 'log');

%% Load psychophysics

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% add additional data from Jonathan, but keep only AC_1_2 plane
% pp_extra = open('data/extra_ternary_thresholds.mat');
% pp_extra_AC12 = selectMeasurements(pp_extra.avg, ...
%     strcmp(pp_extra.avg.groups, 'AC_1_2'));
% 
% pp = catMeasurements(pp, pp_extra_AC12);

if symmetrizePP
    % make sure data is symmetric
    ppOriginal = pp;
    
    reflectTrafo = @(group) applyGroupReflection(group, 3);
    ppReflected = applyToThresholds(pp, reflectTrafo, 'closed', true);
    
    pp = averageMeasurements(ppOriginal, ppReflected);
end

%% Load the multi-cutoff data

ni = open(fullfile('save', 'TernaryDistribution_PennNoSky_2x32_multicutoff_predictions.mat'));

%% Compare to PP

groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;

compType = 'direct';
opts = {'groupMaskFct', groupMaskFct, 'hiLoRatioLimit', 2.0};

nCutoffs = length(ni.allPredictionResults);
comparisons = zeros(nCutoffs, 1);
details = cell(1, nCutoffs);
for i = 1:nCutoffs
    [comparisons(i), details{i}] = compareMeasurements(...
        ni.allPredictionResults(i).predictions, pp, ...
        compType, opts{:});
end

%% Make figure

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [5.08 2];

% ax = axes;
% ax.Position = [0 0 1 1];
% ax.Units = 'inches';
% ax.Position = [0.4 0.5 2.5 0.9];

[colors, colorDict] = get_palette();

hold on;
plot([0, 1], [0 0], 'color', lighten(colorDict('gray'), 0.65), 'linewidth', 0.5);

switch diffType
    case 'spe'
        allDifferencesMatrix = cell2mat(cellfun(@(s) ...
            2*(s.common.measurements2.thresholds - s.common.measurements1.thresholds) ./ ...
              (s.common.measurements2.thresholds + s.common.measurements1.thresholds), ...
            details, 'uniform', false));
    case 'log'
        allDifferencesMatrix = cell2mat(cellfun(@(s) s.common.logdiff, details, ...
            'uniform', false));
    otherwise
        error('Unknown diffType.');
end

allDifferences = allDifferencesMatrix(:);
grayAmounts = arrayfun(@(s) diff(s.cutoffs), ni.allPredictionResults);
allCategories = flatten(meshgrid(grayAmounts(:)', ones(size(allDifferencesMatrix, 1), 1)));
% make sure categories are sorted, to avoid issues
[allCategories, reorder] = sort(allCategories);
grayAmounts = sort(grayAmounts);
allDifferences = allDifferences(reorder);
allColors = repmat(lighten(colorDict('gray'), 0.3), length(grayAmounts), 1);
[~, idxMain] = min(abs(grayAmounts - 1/3));
allColors(idxMain, :) = colorDict('blue');
width = 0.03;
switch plotType
    case 'jitter'
        h = stripPlot(allCategories, allDifferences, 'jitter', width, 'marker', '.', ...
            'sizes', 1e-3, 'kde', true, 'colors', allColors);
    case 'violin'
        h = violinPlot(allCategories, allDifferences, 'width', width, 'colors', allColors);
end
% set(h, 'markerfacealpha', 0.5, 'markeredgealpha', 0.5);

hold on;
% boxplot(allDifferencesMatrix, 'colors', 'k', 'whisker', 0, 'symbol', '');

% boxplot messes up the position!
% ax.Position = [0.4 0.5 2.5 0.9];
ax = gca;

beautifygraph(ax);

% NRstrings = cellfun(@(nr) [int2str(nr(1)) 'x' int2str(nr(2))], valuesNR, 'uniform', false);
% set(ax, 'xtick', 1:length(valuesNR), 'xticklabel', NRstrings, 'xticklabelrotation', 45, ...
%     'yminortick', 'on');
% ax.XAxis.FontSize = 8;
% ax.YAxis.FontSize = 8;

% yl = max(ylim);
yl = 2;
ylim([-yl yl]);

xlabel('Fraction gray');
switch diffType
    case 'spe'
        ylabel('Relative error');
    case 'log'
        ylabel('Log error');
    otherwise
        error('Unknown diffType.');
end

set(ax, 'xminortick', 'off');

preparegraph;

set(ax, 'color','none');
set(fig, 'color', 'none');

safePrint(fullfile('figs', 'draft', 'ternarizationEffects'));

%% SCRATCH

grayAmounts = arrayfun(@(s) diff(s.cutoffs), ni.allPredictionResults);
figure;
plot(grayAmounts, comparisons);
