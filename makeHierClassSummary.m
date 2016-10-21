function makeHierClassSummary(stats, varargin)
% makeHierClassSummary Make a summary figure showcasing the comparison
% between natural and manmade texture statistics in an image database.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% length of confidence interval, as a multiple of standard deviation
parser.addParameter('ciMultiplier', 2, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% figure out a size for the figure
screenPos = get(0,'screensize');
screenW = screenPos(3);
screenH = screenPos(4);

figH = floor(min(screenW, screenH)*0.8);
figW = figH;

figure;
set(gcf, 'units', 'pixels', 'position', ...
    [(screenW-figW)/2 (screenH-figH)/2 figW figH]);
%preparegraph;
papW = figW/72;
papH = figH/72;
set(gcf, 'paperunits', 'inches', 'papersize', [papW papH], ...
    'paperposition', [0 0 papW papH]);
set(gcf, 'PaperPositionMode', 'auto');

% draw the principal component bar plot
subplot('position', [0.05 0.5 0.45 0.4]);
pcAllConds = [stats.summary.both.topPC(:) ...
              stats.summary.natural.topPC(:) ...
              stats.summary.manmade.topPC(:)];
h = bar(pcAllConds, 'grouped');
set(h(1), 'facecolor', [0.5 0.5 0.5]);
set(h(2), 'facecolor', [1 0 0]);
set(h(3), 'facecolor', [1 0.6 0]);
ylim([-1 1]);

% prettify the plot
set(gca, 'box', 'off', 'xaxislocation', 'origin', 'xtick', [], 'tickdir', 'out');
set(gca, 'ytick', [-1 -0.5 0 0.5 1]);
set(gca, 'ticklength', 2*get(gca, 'ticklength'));

% add error bars
hold on;
errPcAllConds = params.ciMultiplier*...
    [stats.summary.both.stdTopPC(:) ...
     stats.summary.natural.stdTopPC(:) ...
     stats.summary.manmade.stdTopPC(:)];
% solution from MathWorks,
% https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
n = size(pcAllConds, 1);
groupWidth = min(0.8, 3/(3+1.5));
for i = 1:3
    x = (1:n) - groupWidth/2 + (2*i-1) * groupWidth / (2*3);
    errorbar(x, pcAllConds(:, i), errPcAllConds(:, i), 'k', ...
        'linestyle', 'none', 'capsize', 3, 'linewidth', 1);
end

text(1.5, 1.2, 'Top principal component', 'fontunits', 'normalized', 'fontsize', 0.07, ...
    'verticalalignment', 'top');

% draw the discriminant bar plot
subplot('position', [0.55 0.5 0.43 0.4]);
discNorm = stats.summary.lda.discriminant/norm(stats.summary.lda.discriminant);
bar(discNorm, 'edgecolor', [1 0.3 0], 'facecolor', 'none', 'linewidth', 1);
ylim([-1 1]);

% prettify the plot
set(gca, 'box', 'off', 'xaxislocation', 'origin', 'xtick', [], 'tickdir', 'out');
set(gca, 'ytick', [-1 -0.5 0 0.5 1]);
set(gca, 'ticklength', 2*get(gca, 'ticklength'));

% add error bars
hold on;
discErr = params.ciMultiplier*stats.summary.lda.stdDiscriminant/norm(stats.summary.lda.discriminant);
errorbar(discNorm, discErr, 'k', 'linestyle', 'none', 'capsize', 3, ...
    'linewidth', 1);

text(0.5, 1.2, 'Discriminant', 'fontunits', 'normalized', 'fontsize', 0.07, ...
    'verticalalignment', 'top');
text(5, 1.2, 'Confusion matrix', 'fontunits', 'normalized', 'fontsize', 0.07, ...
    'verticalalignment', 'top');

% show confusion matrices based on LDA
axes('position', [0.78 0.8 0.07 0.07]);
drawStochasticMatrix(stats.summary.lda.trainingContingency, 'edgeColor', 'w');
text(1.5, -0.1, '(in sample)', 'fontunits', 'normalized', 'fontsize', 0.2, ...
    'horizontalalignment', 'center', 'verticalalignment', 'top');

axes('position', [0.88 0.8 0.07 0.07]);
drawStochasticMatrix(stats.summary.lda.validationContingency, 'edgeColor', 'w');
text(1.5, -0.1, '(out of sample)', 'fontunits', 'normalized', 'fontsize', 0.2, ...
    'horizontalalignment', 'center', 'verticalalignment', 'top');

% show cosine similarity matrix
subplot('position', [0.16 0.05 0.3 0.3]);
drawSymmetricMatrix(stats.summary.similarity.matrix, [0 1], 'abs', true);
cAx = gca;
colormap(cAx, flipud(colormap('gray')));
%insetPos = get(cAx, 'tightinset');

cBar = colorbar('location', 'SouthOutside');
cAxPos = get(cAx, 'position');
cBarPos = get(cBar, 'position');
%set(cBar, 'position', [cAxPos(1)+insetPos(2) cAxPos(2)-0.018 cAxPos(3)-2*insetPos(2) 0.25*cBarPos(4)]);
set(cBar, 'position', [cAxPos(1) cAxPos(2)-0.018 cAxPos(3) 0.25*cBarPos(4)]);
set(cBar, 'ticklength', 0, 'ticks', [0 1]);
xLabH = xlabel(cBar, 'cosine similarity', 'fontunits', 'normalized', 'fontsize', 2.5);
set(xLabH, 'position', get(xLabH, 'Position') + [0 .6 0]);

hold on;
text(0.3, 1, ['natural+' ; 'manmade '], 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06);
text(0.3, 2, 'natural', 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06, 'color', [1 0 0]);
text(0.3, 3, 'manmade', 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06, 'color', [1 0.6 0]);
text(0.3, 4, 'manmade', 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06, 'color', [1 0.6 0]);
text(-0.85, 4, '/', 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06);
text(-0.95, 4, 'natural', 'verticalalignment', 'middle', 'horizontalalignment', 'right', ...
    'fontunits', 'normalized', 'fontsize', 0.06, 'color', [1 0 0]);

% show standard error
axes('position', [0.32 0.32 0.14 0.14]);
tmp = 10^ceil(log10(max(stats.summary.similarity.stdMatrix(:))));
drawSymmetricMatrix(stats.summary.similarity.stdMatrix, [0 tmp], 'abs', true);
cAxStd = gca;
colormap(cAxStd, flipud(colormap('gray')));

cBarStd = colorbar('location', 'SouthOutside');
cAxStdPos = get(cAxStd, 'position');
cBarStdPos = get(cBarStd, 'position');
%set(cBar, 'position', [cAxPos(1)+insetPos(2) cAxPos(2)-0.018 cAxPos(3)-2*insetPos(2) 0.25*cBarPos(4)]);
set(cBarStd, 'position', [cAxStdPos(1) cAxStdPos(2) - 0.01 cAxStdPos(3) 0.75*cBarStdPos(4)]);
set(cBarStd, 'ticklength', 0, 'ticks', [0 tmp]);
xLabStdH = xlabel(cBarStd, 'std', 'fontunits', 'normalized', 'fontsize', 3);
set(xLabStdH, 'position', get(xLabStdH, 'Position') + [0 2 0]);

% show discriminants in frame generated by top 2 PCs
subplot('position', [0.55 0.05 0.43 0.4]);
hold on;
discrProj = stats.reps.lda.discriminant'*stats.summary.both.allPC;
nreps = size(stats.reps.lda.discriminant, 2);
for i = 1:nreps
    threshPt = -stats.reps.lda.const(i)*discrProj(i, :)/norm(discrProj(i, :))^2;
    line([threshPt(1), -discrProj(i, 1)], [threshPt(2), -discrProj(i, 2)], 'color', [1 0 0]);
    line([threshPt(1), discrProj(i, 1)], [threshPt(2), discrProj(i, 2)], 'linestyle', ':', 'color', [1 0.6 0]);
end
plot(-discrProj(:, 1), -discrProj(:, 2), '.', 'color', [1 0 0]);
plot(discrProj(:, 1), discrProj(:, 2), '.', 'color', [1 0.6 0]);

discrSummProj = stats.summary.lda.discriminant'*stats.summary.both.allPC;
discrSummErrProj = stats.summary.lda.stdDiscriminant'*stats.summary.both.allPC;
txt1Pos = -1.2*discrSummProj(:, 1:2) - discrSummErrProj(:, 1:2);
txt2Pos = 1.2*discrSummProj(:, 1:2) + discrSummErrProj(:, 1:2);
text(txt1Pos(1), txt1Pos(2), 'natural', 'color', [1 0 0], ...
    'horizontalalignment', 'center', 'fontunits', 'normalized', ...
    'fontsize', 0.05);
text(txt2Pos(1), txt2Pos(2), 'manmade', 'color', [1 0.6 0], ...
    'horizontalalignment', 'center', 'fontunits', 'normalized', ...
    'fontsize', 0.05);

xlabel('first principal component');
ylabel('second principal component');

xRng = max(abs(min(txt1Pos(1), txt2Pos(1))), abs(max(txt1Pos(1), txt2Pos(1))));
yRng = max(abs(min(txt1Pos(2), txt2Pos(2))), abs(max(txt1Pos(2), txt2Pos(2))));
set(gca, 'box', 'on');

xlim(1.3*[-xRng xRng]);
ylim(1.3*[-yRng yRng]);

end