% make plot showing the dependence of the NI-PP match on the choice of N
% and R parameters

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   dbChoice
%       Choice of database to use. Available options are
%           'PennNoSky' -- Penn image database, filtering out pictures with
%                          lots of sky
%           'vanHateren'    -- van Hateren image database, filtering out
%                              pictures with lots of sky or lots of
%                              human-made objects
%   compressType
%       Choose the way in which image values were compressed in the [0, 1]
%       interval before ternarizing. Options are
%           'equalize' -- histogram equalization
%           'contrast' -- contrast adaptation
%   valuesNR
%       Cell array of pairs [N, R] of downsampling factor and patch size
%       for which to show NI-PP match.
%   symmetrizePP
%       Set to `true` to take the average between each psychophysics
%       measurement and the measurement in the opposite texture direction.
%       This effectively forces the measurements to be centered at the
%       origin.

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('valuesNR', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});
setdefault('symmetrizePP', false);

%% Preprocess options

valuesN = cellfun(@(x) x(1), valuesNR);
valuesR = cellfun(@(x) x(2), valuesNR);

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end

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

%% Load each NI prediction, compare to PP

groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;
% groupMaskFct = @(g) length(g) > 3;
% groupMaskFct = @(g) true;

compType = 'direct';
opts = {'groupMaskFct', groupMaskFct, 'hiLoRatioLimit', 2.0};

comparisons = zeros(size(valuesNR));
details = cell(size(valuesNR));
for i = 1:length(valuesNR)
    crtNR = valuesNR{i};
    NRstr = [int2str(crtNR(1)) 'x' int2str(crtNR(2))];
    crtFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr '.mat'];
    niStructure = open(fullfile('save', crtFileName));
    ni = niStructure.predictions;
    
    [comparisons(i), details{i}] = compareMeasurements(ni, pp, compType, opts{:});
end

%% Make figure

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [3 1.5];

ax = axes;
ax.Units = 'inches';
ax.Position = [0.4 0.5 2.5 0.9];

hold on;
plot([0, length(valuesNR)+1], [0 0], 'color', [0.65 0.65 0.65], 'linewidth', 0.5);

% allDifferencesMatrix = cell2mat(cellfun(@(s) s.common.logdiff, details, ...
%     'uniform', false));
allDifferencesMatrix = cell2mat(cellfun(@(s) ...
    2*(s.common.measurements2.thresholds - s.common.measurements1.thresholds) ./ ...
      (s.common.measurements2.thresholds + s.common.measurements1.thresholds), ...
    details, 'uniform', false));

allDifferences = allDifferencesMatrix(:);
allCategories = flatten(meshgrid(1:length(valuesNR), ones(size(allDifferencesMatrix, 1), 1)));
h = stripPlot(allCategories, allDifferences, 'jitter', 0.5, 'sizes', 1, 'kde', true);
set(h, 'markerfacealpha', 0.5, 'markeredgealpha', 0.5);

hold on;
boxplot(allDifferencesMatrix, 'colors', 'k', 'whisker', 0, 'symbol', '');

% boxplot messes up the position!
ax.Position = [0.4 0.5 2.5 0.9];

beautifygraph(ax, 'fontscale', 0.6667);

NRstrings = cellfun(@(nr) [int2str(nr(1)) 'x' int2str(nr(2))], valuesNR, 'uniform', false);
set(ax, 'xtick', 1:length(valuesNR), 'xticklabel', NRstrings, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;

% yl = max(ylim);
yl = 2;
ylim([-yl yl]);

% ylabel('Log error');
ylabel('Relative error');

set(ax, 'xminortick', 'off');

preparegraph;

safePrint(fullfile('figs', 'draft', 'NRdependence'));

%% Make figure (old)

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 3 1.4];

h = bar(comparisons, 0.6);

h.FaceColor = [0 0.3438 0.7410];
h.EdgeColor = 'none';

ax = gca;
NRstrings = cellfun(@(nr) [int2str(nr(1)) 'x' int2str(nr(2))], valuesNR, 'uniform', false);
set(ax, 'xtick', 1:length(valuesNR), 'xticklabel', NRstrings, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ax.XAxis.FontSize = 8;
ylabel('RMS log error');

beautifygraph('fontscale', 0.6667);
preparegraph;

set(ax, 'xminortick', 'off');

safePrint(fullfile('figs', 'draft', 'NRdependence'));