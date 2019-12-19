% make plots comparing NI predictions to measured PP thresholds for ternary
% textures

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the NI predictions

load(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32_square.mat'));

%% Single planes

rng(133);
plusMinus = '+-';
plotTernaryMatrix({predictions, ternaryAvg}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'edgelabels', 'none', 'axisovershoot', 0}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedShape', [2 2], 'fixedSize', [4 4], ...
    'fixedSizeUnits', 'inches'}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 10, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'drawPatches', true, ...
    'patchOptionsSingle', {'patchExtent', 0.34, 'offset', 0.25});
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternarySecondOrderNIPPMatch.pdf'));

%% Mixed planes

fig = figure;
fig.Units = 'inches';
totalX = 6;
totalY = 6;
fig.Position = [2 2 totalX totalY];

ax = zeros(22, 1);
figX = 1.16;
figY = 1.16;
factorX = 0.8;
factorY = 0.8;
edgeX = 0.10;
edgeY = 0.05;
for i = 1:length(ax)
    crtAx = axes;
    crtRow = floor((i-1)/5);
    crtCol = mod(i-1, 5);
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [crtCol*figX + edgeX totalY - (crtRow+1)*figY - edgeY ...
        figX*factorX figY*factorY];
    ax(i) = crtAx;
end

plusMinus = '+-';
plotTernaryMatrix({predictions, ternaryAvg}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', ax}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.42, ...
        'coeffToStr', @(i) plusMinus(i), 'squareSubscripts', true}, ...
    'xLabelAlignment', {'center', 'bottom'}, 'yLabelAlignment', {'left', 'bottom'});
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternaryMixedNIPPMatch.pdf'));

%% Mixed planes (OLD)

fig = figure;
fig.Units = 'inches';
totalX = 6;
totalY = 5;
fig.Position = [2 2 totalX totalY];

ax = zeros(22, 1);
rowNumbers = [2 2 6 6 6];
rowShifts = [4 4 0 0 0];
rowStarts = [1 1+cumsum(rowNumbers(1:end-1))];
figX = totalX/max(rowNumbers);
figY = totalY/length(rowNumbers);
factorX = 0.75;
factorY = 0.75;
for i = 1:length(ax)
    crtAx = axes;
    crtRow = find(i >= rowStarts, 1, 'last');
    crtShift = rowShifts(crtRow);
    crtCol = i - rowStarts(crtRow); % 0-based!!
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [(crtShift + crtCol)*figX totalY - crtRow*figY figX*factorX figY*factorY];
    ax(i) = crtAx;
end

plusMinus = '+-';
plotTernaryMatrix({predictions, ternaryAvg}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', ax}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'xLabelAlignment', {'center', 'bottom'}, 'yLabelAlignment', {'left', 'bottom'});
preparegraph;
