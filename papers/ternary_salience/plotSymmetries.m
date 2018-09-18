% make plot showing effect of symmetry transformations

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
%   symmetryNR
%       Choose one of the analyses for the symmetry plot, based on
%       block-averaging factor (N) and patch size (R).
%   restrictToFocus
%       Set to `true` to only keep patches that were identified as in-focus
%       by a two-Gaussian fit.
%   symmetrizePP
%       Set to `true` to take the average between each psychophysics
%       measurement and the measurement in the opposite texture direction.
%       This effectively forces the measurements to be centered at the
%       origin.

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('symmetryNR', [2, 32]);
setdefault('restrictToFocus', true);
setdefault('symmetrizePP', false);

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];
NRstr = [int2str(symmetryNR(1)) 'x' int2str(symmetryNR(2))];
symmetryNIPredFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr '.mat'];

%% Load the data

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

%% Load the ternary NI distribution

niStatsAll = open(fullfile('save', niFileName));

% choose one of the analyses
idx = find(cellfun(@(nr) isequal(nr, symmetryNR), niStatsAll.valuesNR));
if isempty(idx)
    error('Can''t find NR selection in NI distribution file.');
end
if length(idx) > 1
    error('Found multiple matches to the NR selection.');
end

niStats0 = niStatsAll.results{idx};

% check that the distribution we loaded has focus information
niStats = rmfield(niStats0, 'focus');
if restrictToFocus && isfield(niStats0, 'focus')
    disp('Restricting to in-focus patches.');
    mask = (niStats0.focus.clusterIds == niStats0.focus.focusCluster);
    fields = {'ev', 'patchLocations', 'imageIds'};
    for i = 1:length(fields)
        niStats.(fields{i}) = niStats.(fields{i})(mask, :);
    end
    niStats.covM = cov(niStats.ev);
end

%% Load the NI predictions

load(fullfile('save', symmetryNIPredFileName));
ni = predictions;

%% Calculate the effects for a series of transformations

% set up the transformations we want to look at 
trafos = containers.Map;
trafos('lrFlip') = @(group) applyGroupGeometricPermutation(group, 3, 'BADC');
trafos('udFlip') = @(group) applyGroupGeometricPermutation(group, 3, 'CDAB');
trafos('rot90') = @(group) applyGroupGeometricPermutation(group, 3, 'CADB');
trafos('rot180') = @(group) applyGroupGeometricPermutation(group, 3, 'DCBA');
trafos('rot270') = @(group) applyGroupGeometricPermutation(group, 3, 'BDAC');
trafos('exch12') = @(group) applyGroupColorTransformation(group, 3, 2, 0);
trafos('exch01') = @(group) applyGroupColorTransformation(group, 3, 2, 1);
trafos('exch02') = @(group) applyGroupColorTransformation(group, 3, 2, 2);
trafos('cycle120') = @(group) applyGroupColorTransformation(group, 3, 1, 1);
trafos('cycle201') = @(group) applyGroupColorTransformation(group, 3, 1, 2);

% generate the transformed NI predictions & PP data
trafoNames = trafos.keys;
nTrafos = length(trafoNames);

transformedNI = cell(1, nTrafos);
transformedPP = cell(1, nTrafos);

for i = 1:nTrafos
    trafoName = trafoNames{i};
    fct = trafos(trafoName);
    
    disp(['Working on ' trafoName '...']);
    
    % apply to NI stats
    niStatsTransformed = niStats;
    niStatsTransformed.ev = applyToStats(niStats.ev, 3, fct);
    
    % apply to PP thresholds
    transformedPP{i} = applyToThresholds(pp, fct, 'closed', true);
    
    % find transformed NI threshold predictions
    [~, transformedNI{i}] = getPredictionsFromTernaryStats(niStatsTransformed.ev, pp, ...
        'fitScaleOptions', {'fixScale', predictionDetails.aCoeff}, ...
        'efficientCodingOptions', {'gainTransform', gainTransformFct});
    
    % find which directions changed and which didn't in the NI
    [~, predictionsShuffle] = applyToThresholds(ni, fct);
    transformedNI{i}.changed = flatten(predictionsShuffle(:)' ~= 1:length(predictionsShuffle));
end

%% Check how large the effect of each transformation was on NI & PP

compType = 'direct';
niEffectSizes = zeros(nTrafos, 1);
ppEffectSizes = zeros(nTrafos, 1);
niChangedCounts = zeros(nTrafos, 1);
ppChangedCounts = zeros(nTrafos, 1);

% groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;
groupMaskFct = @(g) length(g) > 3;
% groupMaskFct = @(g) true;

opts = {'groupMaskFct', groupMaskFct, 'hiLoRatioLimit', 2.0};

for i = 1:nTrafos
    [niEffectSizes(i), niDetails] = compareMeasurements(...
        ni, transformedNI{i}, compType, opts{:});
    [ppEffectSizes(i), ppDetails] = compareMeasurements(...
        pp, transformedPP{i}, compType, opts{:});
    
    if strcmp(compType, 'group') || strcmp(compType, 'ellipse')
        niChangedCounts(i) = numel(niDetails.common.uniqueGroups);
        ppChangedCounts(i) = numel(ppDetails.common.uniqueGroups);
    elseif strcmp(compType, 'direct')
        niChangedCounts(i) = niDetails.common.nAveraged;
        ppChangedCounts(i) = ppDetails.common.nAveraged;
    end
end

%% Make figure

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 3 1.4];

plotMax = max(max(niEffectSizes, max(ppEffectSizes)))*1.2;
h = bar([niEffectSizes(:) ppEffectSizes(:)], 0.9);

h(1).FaceColor = [0 0.3438 0.7410];
h(2).FaceColor = [0.8000 0.3000 0.3000];

h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';

ax = gca;
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ax.XAxis.FontSize = 8;
ylabel('Change');
% ylim([0, 0.45]);
ylim([0, plotMax]);
xlim([0, nTrafos+1]);

legend({'NI', 'PP'}, 'fontsize', 8);

beautifygraph('fontscale', 0.6667);
preparegraph;

set(ax, 'xminortick', 'off');

safePrint(fullfile('figs', 'draft', 'symmetries'));

%% SCRATCH

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 8 6];

plotMax = max(max(niEffectSizes, max(ppEffectSizes)))*1.2;

ax = axes;
ax.OuterPosition = [0 0.5 1 0.5];
% subplot(2, 1, 1);
bar(niEffectSizes);
title('Natural images');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('Average \Delta(log thresholds)');
% ylim([0, 0.45]);
ylim([0, plotMax]);
xlim([0, nTrafos+1]);

for i = 1:size(niChangedCounts, 1)
    text(i, niEffectSizes(i), int2str(niChangedCounts(i)), ...
        'verticalalignment', 'bottom', 'horizontalalignment', 'center');
end

ax = axes;
ax.OuterPosition = [0 0 1 0.5];
% subplot(2, 1, 2);
bar(ppEffectSizes);
title('Psychophysics');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('Average \Delta(log thresholds)');
% ylim([0, 0.45]);
ylim([0, plotMax]);
xlim([0, nTrafos+1]);

for i = 1:size(ppChangedCounts, 1)
    text(i, ppEffectSizes(i), int2str(ppChangedCounts(i)), ...
        'verticalalignment', 'bottom', 'horizontalalignment', 'center');
end
