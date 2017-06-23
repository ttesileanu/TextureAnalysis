% show how texture stats change upon using continuous stats

%% Generate results

filter_file = open('Ffilter1x32.mat');
filter_file_AF2 = open('Ffilter2x32.mat');
filter = filter_file.Ffilter;
filter_AF2 = filter_file_AF2.Ffilter;

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

tic;
res0 = analyzeImageSet(images, 'NaturalImages', 1, filter, 2);
disp(['Binary stats took ' num2str(toc, '%.2f') ' seconds to generate from scratch.']);

tic;
res1 = analyzeImageSet(images, 'NaturalImages', 1, filter, inf);
disp(['Continuous stats took ' num2str(toc, '%.2f') ' seconds to generate from scratch.']);

tic;
res0_AF2 = analyzeImageSet(images, 'NaturalImages', 2, filter_AF2, 2);
disp(['Binary stats, blockAF=2 took ' num2str(toc, '%.2f') ' seconds to generate from scratch.']);

tic;
res1_AF2 = analyzeImageSet(images, 'NaturalImages', 2, filter_AF2, inf);
disp(['Continuous stats, blockAF=2 took ' num2str(toc, '%.2f') ' seconds to generate from scratch.']);

%% Save results

save('binary_vs_continuous.mat', 'res0', 'res1', 'res0_AF2', 'res1_AF2', ...
    'filter', 'filter_AF2', 'images');

%% Load results

load('binary_vs_continuous.mat');

%% Make plots

fig = figure;
fig.Units = 'inches';
fig.Position = [1 2 18 6];

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

for i = 1:10
    subplot(2, 5, i);

    hold on;
    
    plot([0 0], [-1 1], '--', 'color', [0.3 0.3 0.3]);
    plot([-1 1], [0 0], '--', 'color', [0.3 0.3 0.3]);

    smartscatter(res0.ev(:, i), res1.ev(:, i), 'r.', ...
        'alpha', 0.1, 'maxPoints', 50000);
    
    max_ax = max([abs(res0.ev(:, i)) ; abs(res1.ev(:, i))]);
    xlim([-max_ax max_ax]);
    ylim([-max_ax max_ax]);

    axis equal;
    xlabel('binary');
    ylabel('continuous');
    title(labels{i});
    beautifygraph;
    
    drawfitline(res0.ev(:, i), res1.ev(:, i), 'style', {'k', 'linewidth', 1});
end

preparegraph;

%% Make plots for blockAF = 2

fig = figure;
fig.Units = 'inches';
fig.Position = [1 2 18 6];

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

for i = 1:10
    subplot(2, 5, i);

    hold on;
    
    plot([0 0], [-1 1], '--', 'color', [0.3 0.3 0.3]);
    plot([-1 1], [0 0], '--', 'color', [0.3 0.3 0.3]);

    smartscatter(res0_AF2.ev(:, i), res1_AF2.ev(:, i), 'r.', ...
        'alpha', 0.1, 'maxPoints', 50000);
    
    max_ax = max([abs(res0_AF2.ev(:, i)) ; abs(res1_AF2.ev(:, i))]);
    xlim([-max_ax max_ax]);
    ylim([-max_ax max_ax]);

    axis equal;
    xlabel('binary');
    ylabel('continuous');
    title(labels{i});
    beautifygraph;
    
    drawfitline(res0_AF2.ev(:, i), res1_AF2.ev(:, i), 'style', {'k', 'linewidth', 1});
end

preparegraph;