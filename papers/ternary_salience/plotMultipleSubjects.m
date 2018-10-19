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
totalY = 7.5;
fig.Position = [2 2 totalX totalY];

ax = zeros(26, 1);
figX = 1.2;
figY = 1.2;
figXSingle = 1.5;
figYSingle = 1.5;
factorX = 0.8;
factorY = 0.8;
edgeX = 0.05;
edgeY = -0.1;
crtX = edgeX + (figXSingle - figX)/2;
crtY = totalY - figY - edgeY;
for i = 1:length(ax)
    crtAx = axes;
        
    crtAx.Units = 'inches';
%     crtAx.OuterPosition = [crtCol*figX + edgeX totalY - (crtRow+1)*figY - edgeY ...
%         figX*factorX figY*factorY];
    crtAx.OuterPosition = [crtX crtY figX*factorX figY*factorY];
    ax(i) = crtAx;
    
    % keep only single-group planes on first row
    if i < 4
        crtX = crtX + figXSingle;
    else
        if i == 4
            crtX = edgeX;
            crtY = crtY - figYSingle;
        else
            crtX = crtX + figX;
            if crtX > totalX
                crtX = edgeX;
                crtY = crtY - figY;
            end
        end
    end
end

% plot single planes
plusMinus = '+-';
plotTernaryMatrix(pp, 'ellipse', false, ...
    'markers', {'x'}, 'sizes', 2, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', ax(1:4)}, ...
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
    'plotterOptions', {'fixedAxes', ax(5:end)}, ...
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
