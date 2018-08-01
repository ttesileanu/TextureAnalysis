% make plots comparing NI predictions to measured PP thresholds for ternary
% textures

%% Load the data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the NI predictions

load(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32.mat'));

%% Single planes

plotTernaryMatrix({predictions, ternaryAvg}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedShape', [1 4], 'fixedSize', [6 1.5], ...
    'fixedSizeUnits', 'inches'});
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternarySecondOrderNIPPMatch.pdf'));

%% Mixed planes

fig = figure;
fig.Units = 'inches';
fig.Position = [2 2 6 5];

ax = zeros(22, 1);
% row_numbers = [3 3 5 6 5];
% row_shifts = [3 3 0.5 0 0.5];
rowNumbers = [2 2 6 6 6];
rowShifts = [4 4 0 0 0];
rowStarts = [1 1+cumsum(rowNumbers(1:end-1))];
figX = 1/max(rowNumbers);
figY = 1/length(rowNumbers);
for i = 1:length(ax)
    crtAx = axes;
    crtRow = find(i >= rowStarts, 1, 'last');
    crtShift = rowShifts(crtRow);
    crtCol = i - rowStarts(crtRow); % 0-based!!
    crtAx.OuterPosition = [(crtShift + crtCol)*figX 1 - crtRow*figY figX figY];
    ax(i) = crtAx;
end

plotTernaryMatrix({predictions, ternaryAvg}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', ax});
preparegraph;

safePrint(fullfile('figs', 'draft', 'ternaryMixedNIPPMatch.pdf'));