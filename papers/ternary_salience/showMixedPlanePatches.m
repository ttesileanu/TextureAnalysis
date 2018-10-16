% try to draw texture patches for every mixed plane

% make plots comparing NI predictions to measured PP thresholds for ternary
% textures

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Draw the textures

rng(82);

fig = figure;
fig.Units = 'inches';
totalX = 6;
totalY = 8.5;
fig.Position = [2 2 totalX totalY];

mixedGroups = sortGroups(unique(ternaryAvg.groups(cellfun(@(s) sum(s == ';') > 0, ...
    ternaryAvg.groups))));
ax = zeros(length(mixedGroups), 1);
figX = 1.42;
figY = 1.42;
edgeX = -0.1;
edgeY = 0.22;
plusMinus = '+-';
for i = 1:length(ax)
    crtAx = axes;
    crtRow = floor((i-1)/4);
    crtCol = mod(i-1, 4);
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [crtCol*figX + edgeX totalY - (crtRow+1)*figY - edgeY ...
        figX figY];
    ax(i) = crtAx;
    
    hold on;
%     drawTernaryMixedBackground('circles', []);
    showTernaryTexturesMixed(mixedGroups{i}, ...
        'patchExtent', 0.27, 'fixDistance', 0.45, ...
        'saturation', 1.0, ...
        'saturationScale', 'relative', ...
        'drawRadii', true);
    
    beautifygraph;
    
    axis equal;
    axis off;
    
    xlim([-0.6, 0.6]);
    ylim([-0.6, 0.6]);
    
    groupNames = strtrim(split(mixedGroups{i}, ';'));
    textTexGroup(0.59, 0.1, ...
        groupNames{1}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'coeffToStr', @(i) plusMinus(i), ...
        'FontSize', 8);
    textTexGroup(0, 0.59, ...
        groupNames{2}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'coeffToStr', @(i) plusMinus(i), ...
        'FontSize', 8);
end
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternaryMixedPatches.pdf'));
