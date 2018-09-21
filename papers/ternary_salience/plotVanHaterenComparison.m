% make NI-PP match plot for van Hateren database

%% Load the data

[pp, ppRaw] = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the NI predictions

load(fullfile('save', 'TernaryNIPredictions_vanHateren_2x32_square.mat'));

%% Plot the data

uniqueSubjects = unique(pp.subjects);
uniqueColors = lines(length(uniqueSubjects));
colorMap = containers.Map(uniqueSubjects, num2cell(uniqueColors, 2));

fig = figure;
fig.Units = 'inches';
totalX = 6;
totalY = 5;
fig.Position = [2 2 totalX totalY];

% axes for mixed planes
multi_ax = zeros(22, 1);
rowNumbers = [2 2 6 6 6];
rowShifts = [4 4 0 0 0];
rowStarts = [1 1+cumsum(rowNumbers(1:end-1))];
figX = totalX/max(rowNumbers);
figY = totalY/length(rowNumbers);
factorX = 0.75;
factorY = 0.75;
for i = 1:length(multi_ax)
    crtAx = axes;
    
    crtRow = find(i >= rowStarts, 1, 'last');
    crtShift = rowShifts(crtRow);
    crtCol = i - rowStarts(crtRow); % 0-based!!
    
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [(crtShift + crtCol)*figX totalY - crtRow*figY figX*factorX figY*factorY];
    
    multi_ax(i) = crtAx;
end

% axes for single planes
single_ax = zeros(4, 1);
for i = 1:length(single_ax)
    crtAx = axes;
    
    crtRow = floor((i-1)/2);
    crtCol = mod(i-1, 2);
    
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [(2/3 + 5/3*crtCol)*figX totalY - (1 + crtRow)*figY figX figY];
    
    single_ax(i) = crtAx;
end

% plot single planes
plusMinus = '+-';
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', single_ax}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)});

% plot mixed planes
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', multi_ax}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'xLabelAlignment', {'center', 'bottom'}, ...
    'yLabelAlignment', {'left', 'bottom'});
preparegraph;

set(fig, 'Renderer', 'painters')

safePrint(fullfile('figs', 'draft', 'vanHaterenMatch'));
