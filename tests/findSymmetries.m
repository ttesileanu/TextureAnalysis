% find symmetries in the NI predictions and in the PP data

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   dbChoice
%       Choice of database to use. Available options are
%           'PennNoSky' -- Penn image database, filtering out pictures with
%                          lots of sky
%   compressType
%       Choose the way in which image values were compressed in the [0, 1]
%       interval before ternarizing. Options are
%           'equalize' -- histogram equalization
%           'contrast' -- contrast adaptation
%   NRselection
%       Choose one of the analyses, based on block-averaging factor (N) and
%       patch size (R).
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
setdefault('NRselection', [2, 32]);
setdefault('restrictToFocus', true);
setdefault('symmetrizePP', false);

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];
NRstr = [int2str(NRselection(1)) 'x' int2str(NRselection(2))];
niPredFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr '.mat'];

%% Load the data

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% add additional data from Jonathan, but keep only AC_1_2 plane
pp_extra = open('data/extra_ternary_thresholds.mat');
pp_extra_AC12 = selectMeasurements(pp_extra.avg, ...
    strcmp(pp_extra.avg.groups, 'AC_1_2'));

pp = catMeasurements(pp, pp_extra_AC12);

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
idx = find(cellfun(@(nr) isequal(nr, NRselection), niStatsAll.valuesNR));
if isempty(idx)
    error('Can''t find NR selection in NI distribution file.');
end
if length(idx) > 1
    error('Found multiple matches to the NR selection.');
end

niStats0 = niStatsAll.results{idx};

%% NI distribution preprocessing

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

load(fullfile('save', niPredFileName));
ni = predictions;

%% Calculate the effects for a series of transformations

% set up the transformations we want to look at 
trafos = containers.Map;
trafos('lrFlip') = @(group) applyGroupGeometricPermutation(group, 3, 'BADC');
trafos('udFlip') = @(group) applyGroupGeometricPermutation(group, 3, 'CDAB');
trafos('rot90') = @(group) applyGroupGeometricPermutation(group, 3, 'CADB');
trafos('rot180') = @(group) applyGroupGeometricPermutation(group, 3, 'DCBA');
trafos('rot270') = @(group) applyGroupGeometricPermutation(group, 3, 'BDAC');
trafos('exchange12') = @(group) applyGroupColorTransformation(group, 3, 2, 0);
trafos('exchange01') = @(group) applyGroupColorTransformation(group, 3, 2, 1);
trafos('exchange02') = @(group) applyGroupColorTransformation(group, 3, 2, 2);
trafos('cycle120') = @(group) applyGroupColorTransformation(group, 3, 1, 1);
trafos('cycle201') = @(group) applyGroupColorTransformation(group, 3, 1, 2);
trafos('reflect') = @(group) applyGroupReflection(group, 3);
trafos('lr__x02') = @(group) chainGroupTransformations(...
    trafos('lrFlip'), trafos('exchange02'), group);
trafos('ud__x02') = @(group) chainGroupTransformations(...
    trafos('udFlip'), trafos('exchange02'), group);
trafos('r90__x02') = @(group) chainGroupTransformations(...
    trafos('rot90'), trafos('exchange02'), group);
trafos('r180__x02') = @(group) chainGroupTransformations(...
    trafos('rot180'), trafos('exchange02'), group);
trafos('r270__x02') = @(group) chainGroupTransformations(...
    trafos('rot270'), trafos('exchange02'), group);
trafos('x01__refl') = @(group) chainGroupTransformations(...
    trafos('reflect'), trafos('exchange01'), group);
trafos('x01__reflAB11') = @(group) chainGroupTransformations(...
    @(group) applyGroupReflection(group, 3, 'restrict', @(g) strcmp(g, 'AB_1_1')), ...
    trafos('exchange01'), group);
trafos('x01__reflAB11_12') = @(group) chainGroupTransformations(...
    @(group) applyGroupReflection(group, 3, 'restrict', @(g) ismember(g, {'AB_1_1', 'AB_1_2'})), ...
    trafos('exchange01'), group);
trafos('x01__reflAB12') = @(group) chainGroupTransformations(...
    @(group) applyGroupReflection(group, 3, 'restrict', @(g) strcmp(g, 'AB_1_2')), ...
    trafos('exchange01'), group);

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
        'fitScaleOptions', {'fixScale', predictionDetails.aCoeff});
    
    % find which directions changed and which didn't in the NI
    [~, predictionsShuffle] = applyToThresholds(ni, fct);
    transformedNI{i}.changed = flatten(predictionsShuffle(:)' ~= 1:length(predictionsShuffle));
end

%% Check how large the effect of each transformation was on NI & PP

compTypes = {'direct', 'group', 'ellipse'};
nCompTypes = length(compTypes);
niEffectSizes = zeros(nTrafos, nCompTypes);
ppEffectSizes = zeros(nTrafos, nCompTypes);
niChangedCounts = zeros(nTrafos, nCompTypes);
ppChangedCounts = zeros(nTrafos, nCompTypes);

groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;
% groupMaskFct = @(g) length(g) > 3;
% groupMaskFct = @(g) true;

opts = {'groupMaskFct', groupMaskFct, 'hiLoRatioLimit', 2.0};

for i = 1:nTrafos
    for j = 1:nCompTypes
        [niEffectSizes(i, j), niDetails] = compareMeasurements(...
            ni, transformedNI{i}, compTypes{j}, opts{:});
        [ppEffectSizes(i, j), ppDetails] = compareMeasurements(...
            pp, transformedPP{i}, compTypes{j}, opts{:});
        
        if strcmp(compTypes{j}, 'group') || strcmp(compTypes{j}, 'ellipse')
            niChangedCounts(i, j) = numel(niDetails.common.uniqueGroups);
            ppChangedCounts(i, j) = numel(ppDetails.common.uniqueGroups);
        elseif strcmp(compTypes{j}, 'direct')
            niChangedCounts(i, j) = niDetails.common.nAveraged;
            ppChangedCounts(i, j) = ppDetails.common.nAveraged;
        end
    end
end

niEffectSizesCell = num2cell(niEffectSizes', 1);
ppEffectSizesCell = num2cell(ppEffectSizes', 1);
niEffectSizesTable = table(niEffectSizesCell{:}, ...
    'VariableNames', trafoNames, 'RowNames', compTypes);
ppEffectSizesTable = table(ppEffectSizesCell{:}, ...
    'VariableNames', trafoNames, 'RowNames', compTypes);

plotChoice = 'ellipse';
plotChoiceIdx = find(strcmp(compTypes, plotChoice), 1);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 8 6];

plotMax = max(max(niEffectSizes(:, plotChoiceIdx), ...
    max(ppEffectSizes(:, plotChoiceIdx))))*1.2;

ax = axes;
ax.OuterPosition = [0 0.5 1 0.5];
% subplot(2, 1, 1);
bar(niEffectSizes(:, plotChoiceIdx));
title('Natural images');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('\Delta thresholds');
% ylim([0, 0.45]);
ylim([0, plotMax]);
xlim([0, nTrafos+1]);

for i = 1:size(niChangedCounts, 1)
    text(i, niEffectSizes(i, plotChoiceIdx), int2str(niChangedCounts(i, plotChoiceIdx)), ...
        'verticalalignment', 'bottom', 'horizontalalignment', 'center');
end

ax = axes;
ax.OuterPosition = [0 0 1 0.5];
% subplot(2, 1, 2);
bar(ppEffectSizes(:, plotChoiceIdx));
title('Psychophysics');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('\Delta thresholds');
% ylim([0, 0.45]);
ylim([0, plotMax]);
xlim([0, nTrafos+1]);

for i = 1:size(ppChangedCounts, 1)
    text(i, ppEffectSizes(i, plotChoiceIdx), int2str(ppChangedCounts(i, plotChoiceIdx)), ...
        'verticalalignment', 'bottom', 'horizontalalignment', 'center');
end

%% SCRATCH

for i = 1:length(pp.groups)
    groupParts = strtrim(strsplit(pp.groups{i}, ';'));
    sortedParts = sortGroups(groupParts);
    regroup = buildGroupName(sortedParts{:});
    oldGroup = buildGroupName(groupParts{:});
    if ~strcmp(oldGroup, regroup)
        warning('Difference: %s ~= %s', oldGroup, regroup);
    end
end
