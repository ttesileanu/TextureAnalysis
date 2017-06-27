% perform a multiscale analysis

%% Select blockAF and patch sizes

%originalR = 60;
%N_values = 1:6;
originalR = 128;
N_values = [1, 2, 4, 8, 16];

%% Generate multiscale filters

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

filters = cell(1, length(N_values)); %#ok<*NASGU>
for i = 1:length(N_values)
    crtN = N_values(i);
    filters{i} = generateFourierWhitenFilter(images, 'NaturalImages', crtN, round(originalR/crtN));
end

%% Save multiscale filters

for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = round(originalR/crtN);
    filter = filters{i};
    save(fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']), 'filter');
end

%% Load multiscale filters

filters = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = round(originalR/crtN);
    filterFilename = fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']);
    crtFilter = open(filterFilename);
    fields = fieldnames(crtFilter);
    isfilter = @(f) isnumeric(f) && ismatrix(f) && all(size(f) == [crtR crtR]);
    valid = cellfun(@(s) isfilter(crtFilter.(s)), fields);
    if sum(valid) == 0
        error(['Can''t find valid filter data in ' filterFilename '.']);
    elseif sum(valid) > 1
        error(['Don''t know which field to use from ' filterFilename '.']);
    end
    filterField = fields{valid};
    filters{i} = crtFilter.(filterField);
end

%% Generate the multiscale stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    res{i} = analyzeImageSet(images, 'NaturalImages', crtN, filters{i}, inf);
end

%% Save the multiscale stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
save(fullfile('save', ['natural_nosky_multiscale_' int2str(originalR) ...
    '_AF' filtersToString '.mat']), 'res', 'originalR', 'N_values');

%% Load the multiscale stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
load(fullfile('save', ['natural_nosky_multiscale_' int2str(originalR) ...
    '_AF' filtersToString '.mat']));

%% Run cross-scale PCA

% either including all scales, or only AF=2 and up
all_ev = [];
all_ev_ge2 = [];
for i = 1:length(res)
    all_ev = [all_ev res{i}.ev]; %#ok<AGROW>
    if N_values(i) > 1
        all_ev_ge2 = [all_ev_ge2 res{i}.ev]; %#ok<AGROW>
    end
end

[pcs, pc_proj, pc_vars, ~, pc_explained] = pca(all_ev);
[pcs_ge2, pc_proj_ge2, pc_vars_ge2, ~, pc_explained_ge2] = pca(all_ev_ge2);

%% Run per-scale PCA

each_pc = cell(1, length(res));
each_pc_proj = cell(1, length(res));
each_pc_var = cell(1, length(res));
each_pc_explained = cell(1, length(res));
for i = 1:length(res)
    [each_pc{i}, each_pc_proj{i}, each_pc_var{i}, ~, each_pc_explained{i}] = ...
        pca(res{i}.ev);
end

%% Make Scree plots

fig = figure;
fig.Units = 'inches';
fig.Position = [4 0.5 8 8];

% draw cross-scale plot, all AF
ax_cross = axes;
ax_cross.Units = 'normalized';
ax_cross.OuterPosition = [0 0.6 0.5 0.4];

plot(pc_explained, '.-');
hold on;
plot(cumsum(pc_explained), '.-');
grid on;

ints2string0 = @(v) cell2mat(arrayfun(@(n) [',' int2str(n)], v, 'uniform', false));
drop_first = @(s) s(2:end);
ints2string = @(v) ['[' drop_first(ints2string0(v)) ']'];
title(['Cross-scale PCA, N=' ints2string(N_values) ', N*R=' int2str(originalR)]);
legend({'per PC', 'cumulative'});

xlabel('PC#');
ylabel('Variance explained');

ylim([0 100]);

beautifygraph;

% draw cross-scale plot excluding AF1
ax_cross_ge2 = axes;
ax_cross_ge2.Units = 'normalized';
ax_cross_ge2.OuterPosition = [0.5 0.6 0.5 0.4];

plot(pc_explained_ge2, '.-');
hold on;
plot(cumsum(pc_explained_ge2), '.-');
grid on;

title(['Cross-scale PCA, N=' ints2string(N_values(N_values > 1)) ', N*R=' int2str(originalR)]);
legend({'per PC', 'cumulative'});

xlabel('PC#');
ylabel('Variance explained');

ylim([0 100]);

beautifygraph;

% draw within-scale PCA plots
npx = ceil(sqrt(length(res)));
npy = ceil(length(res) / npx);
ax_separated = cell(length(res), 1);
for i = 1:length(res)
    crt_x = mod(i-1, npx);
    crt_y = floor((i-1) / npx);
    
    ax_separated{i} = axes;
    ax_separated{i}.Units = 'normalized';
    % crt_y == 0     --> ax_cross.OuterPosition(2)*(1 - 1/npy)
    % crt_y == npy-1 --> 0
    ax_separated{i}.OuterPosition = [...
        crt_x/npx ax_cross.OuterPosition(2)*(1 - 1/npy)*(1 - crt_y/(npy-1)) ...
        1/npx ax_cross.OuterPosition(2)/npy];

%    subplot(npy, npx, i);
    plot(each_pc_explained{i}, '.-');
    
    hold on;
    plot(cumsum(each_pc_explained{i}), '.-');
    grid on;
    
    title(['Scree plot N=' int2str(N_values(i)) ', R=' int2str(round(originalR/N_values(i)))]);
    legend({'per PC', 'cumulative'}, 'location', 'east');
    
    xlabel('PC#');
    ylabel('Variance explained');
    
    ylim([0 100]);
    
    beautifygraph;
end

preparegraph;

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
print('-dpdf', fullfile('figs', ['pca_scree_' int2str(originalR) '_AF' filtersToString '.pdf']));

%% Make PC plot, cross-scale all AF

makeMultiscalePCAPlot(pcs, pc_explained, N_values);
filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
print('-dpdf', fullfile('figs', ['pca_crossscale_' int2str(originalR) '_AF' filtersToString '.pdf']));

%% Make PC plot,  cross-scale all except AF=1

makeMultiscalePCAPlot(pcs_ge2, pc_explained_ge2, N_values(2:end));
filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values(N_values > 1), 'uniform', false));
print('-dpdf', fullfile('figs', ['pca_crossscale_' int2str(originalR) '_AF' filtersToString '.pdf']));

%% Make PC plot, per scale

makeMultiscalePCAPlot(each_pc, each_pc_explained, N_values);
filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
print('-dpdf', fullfile('figs', ['pca_perscale_' int2str(originalR) '_AF' filtersToString '.pdf']));

%% Diagnostics: how do natural stats compare to random patches?

% generate a number of random images that ends up approximating the number
% of patches we've had from the natural set
n_patches = size(res{1}.ev, 1);
patches_per_row = 8;
patches_per_col = 8;
patches_per_image = patches_per_row*patches_per_col;
n_random_images = round(n_patches / patches_per_image);

random_image_width = originalR * patches_per_row;
random_image_height = originalR * patches_per_col;

res_random = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    res_random{i} = analyzeImageSet(n_random_images, ...
        @(k) rand(random_image_height, random_image_width), ...
        crtN, [], inf, 'patchSize', round(originalR / crtN));
end

%% Save random multiscale stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
save(fullfile('save', ['random_multiscale_' int2str(originalR) ...
    '_AF' filtersToString '.mat']), 'res_random', 'originalR', 'N_values');

%% Load random multiscale stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
load(fullfile('save', ['random_multiscale_' int2str(originalR) ...
    '_AF' filtersToString '.mat']));

%% Diagnostics plots

labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};
plot_pairs = {[1, 10], [2, 3], [4, 5], [6, 7], [8, 9], [2, 10]};

for i = 1:length(res)
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [1 + 0.2*i, 1 + 0.2*i, 5, 3];
    fig.Name = ['N=' int2str(N_values(i)) ', R=' int2str(round(originalR / N_values(i)))];
    
    for j = 1:length(plot_pairs)
        plt_x = mod(j-1, 3)/3;
        plt_y = floor((j-1)/3)/2;
        ax = axes;
        ax.Units = 'normalized';
        ax.OuterPosition = [plt_x 1-plt_y-0.5 1/3 1/2];

        idx1 = plot_pairs{j}(1);
        idx2 = plot_pairs{j}(2);
        
        hold on;
        
        plot([0 0], [-1 1], '--', 'color', [0.3 0.3 0.3]);
        plot([-1 1], [0 0], '--', 'color', [0.3 0.3 0.3]);
        
        hnat = smartscatter(res{i}.ev(:, idx1), res{i}.ev(:, idx2), 10, [0.8 0 0], 'filled', 'alpha', 0.1);
        hrnd = smartscatter(res_random{i}.ev(:, idx1), res_random{i}.ev(:, idx2), 10, [0.6 0.6 0.6], 'filled', 'alpha', 0.1);
        
        xlabel(labels{idx1});
        ylabel(labels{idx2});
        set(get(gca, 'ylabel'), 'rotation', 0);
        
        if max(abs(res{i}.ev(:, idx1))) > 0.5 || max(abs(res{i}.ev(:, idx2))) > 0.5 || ...
                max(abs(res_random{i}.ev(:, idx1))) > 0.5 || max(abs(res_random{i}.ev(:, idx2))) > 0.5
            warning('Plots don''t show full range of variation.');
        end
        
        xlim([-0.5 0.5]);
        ylim([-0.5 0.5]);
        
        beautifygraph;
        
%         if j == 1
%             leg = legend([hnat.hscatter hrnd.hscatter], {'natural', 'random'});
%             leg.Position = [0 0.5 0.2 0.2];
%         end
    end
    preparegraph;
    print('-dpng', '-r300', fullfile('figs', ['multiscale_diagnostics_' ...
        int2str(N_values(i)) 'x' int2str(round(originalR/N_values(i))) '.png']));
end

%% Look at filters

fig = figure;
fig.Units = 'inches';
fig.Position = [1, 1, 5, 3];

for i = 1:length(filters)
    subplot(2, 3, i);
    
    sz2 = round(size(filters{i}, 1)/2);
    rs_filter = circshift(ifft2(filters{i}), [-sz2 -sz2]);
    flt_range = max(abs(rs_filter(:)));
    imagesc(rs_filter, [-flt_range flt_range]);
    
    colormap(redblue);
    axis image;
    
    title(['N=' int2str(N_values(i))]);
    
    beautifygraph;
end
preparegraph;

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
print('-dpng', '-r300', fullfile('figs', ['multiscale_diagnostics_filters_' ...
    int2str(originalR) '_AF' filtersToString '.png']));

%% Compare statistics across scales for the same patch

fig = figure;
fig.Units = 'inches';
n_plots = length(N_values) - 1;
fig.Position = [1 2 18 1.8*n_plots];

for i = 2:length(N_values)
    plt_y = (i-1)/n_plots;
    for j = 1:10
        plt_x = (j-1)/10;
        
        ax = axes;
        ax.Units = 'normalized';
        ax.OuterPosition = [plt_x 1-plt_y 1/10 1/n_plots];

        hold on;
        
        plot([0 0], [-1 1], '--', 'color', [0.3 0.3 0.3]);
        plot([-1 1], [0 0], '--', 'color', [0.3 0.3 0.3]);
        
        res_x = res{1}.ev(:, j);
        res_y = res{i}.ev(:, j);
        smartscatter(res_x, res_y, 10, [0.8 0 0], 'filled', 'alpha', 0.1);
        title(labels{j});
        
        xlabel('N=1');
        ylabel(['N=' int2str(N_values(i))]);
        
        crt_lim = max(max(abs(res_x)), max(abs(res_y)));
        xlim([-crt_lim, crt_lim]);
        ylim([-crt_lim, crt_lim]);
        
        beautifygraph;
    end
end

preparegraph;
print('-dpng', '-r300', fullfile('figs', ['multiscale_interscale_comparison_' ...
    int2str(originalR) '_AF' filtersToString '.png']));