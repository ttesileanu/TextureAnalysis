% make PP plots for multiple subjects

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   showErrorBars
%       Set to `false` to not draw error bars.

setdefault('showErrorBars', true);

%% Load the data

[pp, ppRaw] = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
    'subjects', '*');

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
plusMinus = '+?';
plotTernaryMatrix(pp, 'ellipse', false, ...
    'markers', {'x'}, 'sizes', 2, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', single_ax}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'colorFct', colorMap);

% plot mixed planes
plotTernaryMatrix(pp, 'ellipse', false, ...
    'markers', {'x'}, 'sizes', 2, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', multi_ax}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'xLabelAlignment', {'center', 'bottom'}, ...
    'yLabelAlignment', {'left', 'bottom'}, ...
    'colorFct', colorMap);
preparegraph;

set(fig, 'Renderer', 'painters')

safePrint(fullfile('figs', 'draft', 'thresholdsPerSubject'));

%% Plot legend

fig1 = figure;
fig1.Units = 'inches';
fig1.Position(3:4) = [1 1.2];
hold on;
h = zeros(length(uniqueSubjects), 1);
for i = 1:length(uniqueSubjects)
    h(i) = plot(nan, nan, 'x', 'color', uniqueColors(i, :), 'markersize', 2);
end
legend(h, uniqueSubjects, 'fontsize', 8);

beautifygraph;
axis off;

preparegraph;

safePrint(fullfile('figs', 'draft', 'subjectLegend'));
