% try to draw texture patches for every mixed plane

% make plots comparing NI predictions to measured PP thresholds for ternary
% textures

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Draw the textures

rng(82);

fig = figure;
fig.Units = 'inches';
totalX = 5.49;
totalY = 6.3;
fig.Position = [2 2 totalX totalY];

mixedGroups = sortGroups(unique(ternaryAvg.groups(cellfun(@(s) sum(s == ';') > 0, ...
    ternaryAvg.groups))));
ax = zeros(length(mixedGroups), 1);
figX = 1.32;
figY = 1.00;
edgeX = -0.1;
edgeY = 0.20;
spacingX = 0;
spacingY = 0.05;
plusMinus = '+-';
for i = 1:length(ax)
    crtAx = axes;
    crtRow = floor((i-1)/4);
    crtCol = mod(i-1, 4);
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [
        crtCol*(figX + spacingX) + edgeX ...
        totalY - (crtRow+1)*(figY + spacingY) - edgeY ...
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
    crtAx.Clipping = 'off';
    plot([0 0], [0.65 0.95], 'k', 'linewidth', 0.75);
    plot([-0.05 0 0.05], [0.87 0.95 0.87], 'k', 'linewidth', 0.75);
%     plot(0, 0.95, '^k');
    
    plot([0.65 1.15], [0 0], 'k', 'linewidth', 0.75);
    plot([1.07 1.15 1.07], [-0.05 0 0.05], 'k', 'linewidth', 0.75);
%     plot(1.05, 0, '>k');
    
    textTexGroup(0.61, 0.03, ...
        groupNames{1}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'coeffToStr', @(i) plusMinus(i), ...
        'FontSize', 8, 'squareSubscripts', true);
    textTexGroup(0.06, 0.59, ...
        groupNames{2}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'coeffToStr', @(i) plusMinus(i), ...
        'FontSize', 8, 'squareSubscripts', true);
    
    drawnow;
end
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternaryMixedPatches.pdf'));
