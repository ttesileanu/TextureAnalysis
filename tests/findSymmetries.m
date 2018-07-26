% find symmetries in the NI predictions and in the PP data

%% Load the data

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the ternary NI distribution

niStatsAll = open(fullfile('save', 'TernaryDistribution_PennNoSky.mat'));

% choose one of the analyses
NRselection = [2, 32];
idx = find(cellfun(@(nr) isequal(nr, NRselection), niStatsAll.valuesNR));
if isempty(idx)
    error('Can''t find NR selection in NI distribution file.');
end
if length(idx) > 1
    error('Found multiple matches to the NR selection.');
end

niStats0 = niStatsAll.results{idx};

%% NI distribution preprocessing

% we typically restrict our analysis to patches identified as in-focus
restrictToFocus = true;

% check that the distribution we loaded has focus information
niStats = rmfield(niStats0, 'focus');
if restrictToFocus && isfield(niStats, 'focus')
    disp('Restricting to in-focus patches.');
    mask = (niStats0.focus.clusterIds == niStats0.focus.focusCluster);
    fields = {'ev', 'patchLocations', 'imageIds'};
    for i = 1:length(fields)
        niStats.(fields{i}) = niStats.(fields{i})(mask, :);
    end
    niStats.covM = cov(niStats.ev);
end

%% Load the NI predictions

load(fullfile('save', ['TernaryNIPredictions_PennNoSky_' int2str(NRselection(1)) ...
    'x' int2str(NRselection(2)) '.mat']));
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

compTypes = {'direct', 'group'};
nCompTypes = length(compTypes);
niEffectSizes = zeros(nTrafos, nCompTypes);
ppEffectSizes = zeros(nTrafos, nCompTypes);

groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;

for i = 1:nTrafos
    for j = 1:nCompTypes
        niEffectSizes(i, j) = compareMeasurements(...
            ni, transformedNI{i}, compTypes{j}, 'groupMaskFct', groupMaskFct);
        ppEffectSizes(i, j) = compareMeasurements(...
            pp, transformedPP{i}, compTypes{j}, 'groupMaskFct', groupMaskFct);
    end
end

niEffectSizesCell = num2cell(niEffectSizes', 1);
ppEffectSizesCell = num2cell(ppEffectSizes', 1);
niEffectSizesTable = table(niEffectSizesCell{:}, ...
    'VariableNames', trafoNames, 'RowNames', compTypes);
ppEffectSizesTable = table(ppEffectSizesCell{:}, ...
    'VariableNames', trafoNames, 'RowNames', compTypes);

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 6 4];

ax = axes;
ax.OuterPosition = [0 0.5 1 0.5];
% subplot(2, 1, 1);
bar(niEffectSizes(:, 2));
title('Natural images');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('\Delta thresholds');
ylim([0, 0.4]);

ax = axes;
ax.OuterPosition = [0 0 1 0.5];
% subplot(2, 1, 2);
bar(ppEffectSizes(:, 2));
title('Psychophysics');
set(ax, 'xtick', 1:nTrafos, 'xticklabel', trafoNames, 'xticklabelrotation', 45, ...
    'yminortick', 'on');
ylabel('\Delta thresholds');
ylim([0, 0.4]);

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