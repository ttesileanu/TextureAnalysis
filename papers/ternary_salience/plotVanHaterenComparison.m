% make NI-PP match plot for van Hateren database

%% Load the data

[pp, ppRaw] = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

%% Load the NI predictions

load(fullfile('save', 'TernaryNIPredictions_vanHateren_2x32_square.mat'));

%% Plot the data

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
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', ax(1:4)}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)});

% plot mixed planes
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', ax(5:end)}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'xLabelAlignment', {'center', 'bottom'}, ...
    'yLabelAlignment', {'left', 'bottom'});
preparegraph;

set(fig, 'Renderer', 'painters')

safePrint(fullfile('figs', 'draft', 'vanHaterenMatch'));
