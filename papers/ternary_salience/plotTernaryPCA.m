% look at principal components

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

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('NRselection', [2, 32]);
setdefault('restrictToFocus', true);

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];

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

%% Do PCA

[pcVectors, pcProjection, ~, ~, pcExplained] = pca(...
    expandTextureStats(niStats.ev, 3));
% [pcVectors, pcProjection, ~, ~, pcExplained] = pca(niStats.ev);

%% Scree plot

figure;
plot(pcExplained, 'k.-');

hold on;
plot(cumsum(pcExplained), 'r.-');

xlabel('Principal component');
ylabel('% explained');

ylim([0 100]);

beautifygraph;
preparegraph;

%% Look at individual principal components

figure;
bar(pcVectors(:, 1));

coords = getCoordinateMapping(3, 'analyzeFull');

ticks = [];
labels = {};
lastShape = '';
lastIdx = 0;
for i = 1:length(coords)+1
    if i < length(coords)
        shape = split(coords{i}, '_');
        shape = shape{1};
    else
        shape = '';
    end
    if ~strcmp(shape, lastShape)
        if lastIdx > 0
            ticks = [ticks (i - 1 + lastIdx)/2 - 1/2]; %#ok<AGROW>
            switch lastShape
                case 'A'
                    crtLabel = '\gamma';
                case 'AB'
                    crtLabel = '\beta_{--}';
                case 'AC'
                    crtLabel = '\beta_|';
                case 'BC'
                    crtLabel = '\beta_/';
                case 'AD'
                    crtLabel = '\beta_\\';
                case 'ABC'
                    crtLabel = '\theta_\lceil';
                case 'ABD'
                    crtLabel = '\theta_\rceil';
                case 'ACD'
                    crtLabel = '\theta_\lfloor';
                case 'BCD'
                    crtLabel = '\theta_\rfloor';
                case 'ABCD'
                    crtLabel = '\alpha';
                otherwise
                    error('This shouldn''t happen!');
            end
            labels = [labels {crtLabel}]; %#ok<AGROW>
        end
        lastIdx = i;
        lastShape = shape;
    end
end

set(gca, 'xtick', ticks, 'xticklabels', labels);

beautifygraph('minorticks', 'off');
preparegraph;

%% Check how many PCs have thresholds below 1

load(fullfile('save', ...
    ['TernaryNIPredictions_' dbChoice '_' int2str(NRselection(1)) ...
    'x' int2str(NRselection(2)) '_square.mat']));

PCThresholds = gainsToThresholds(gain, num2cell(pcVectors, 1));

figure;
plot(PCThresholds(1:65));

%% SCRATCH

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [6 4];

allGroups = getCoordinateMapping(3);
plusMinus = '+?';
for i = 1:length(allGroups)
    textTexGroup(i/34, 0.5, allGroups{i}, 'fontsize', 6, ...
        'subscriptSpacing', -0.65, 'coeffToStr', @(i) plusMinus(i));
end
