% make plots comparing NI predictions to measured PP thresholds for third
% and fourth order ternary texture planes

%% Load the data

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
    'subjects', '*');

%% Load the NI predictions

ni = open(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32_square.mat'));

%%

uniqueSubjects = unique(pp.subjects);
uniqueColors = lines(length(uniqueSubjects));
colorMap = containers.Map(uniqueSubjects, num2cell(uniqueColors, 2));

% color for NI
colorMap('subject') = [0 0.3438 0.7410];

fig = figure;
fig.Units = 'inches';
totalX = 6;
totalY = 3;
fig.Position = [2 2 totalX totalY];

ax = zeros(7, 1);
figX = 1.5;
figY = 1.5;
factorX = 0.8;
factorY = 0.8;
edgeX = 0.05;
edgeY = -0.1;
crtX = edgeX;
crtY = totalY - figY - edgeY;
for i = 1:length(ax)
    crtAx = axes;
        
    crtAx.Units = 'inches';
%     crtAx.OuterPosition = [crtCol*figX + edgeX totalY - (crtRow+1)*figY - edgeY ...
%         figX*factorX figY*factorY];
    crtAx.OuterPosition = [crtX crtY figX*factorX figY*factorY];
    ax(i) = crtAx;
    
    if mod(i, 4) == 0
        crtX = edgeX;
        crtY = crtY - figY;
    else
        crtX = crtX + figX;
        if crtX > totalX
            crtX = edgeX;
            crtY = crtY - figY;
        end
    end
end

plusMinus = '+-';
plotTernaryMatrix({ni.predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) > 6 && sum(g == ';') == 0, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', ax}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'colorFct', colorMap);

preparegraph;

set(fig, 'Renderer', 'painters');

safePrint(fullfile('figs', 'draft', 'higherOrderThresholds'));
