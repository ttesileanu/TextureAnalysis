% perform a multiscale analysis

%% Select blockAF and patch sizes

% [N, R] pairs
NR_values = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]};

N_values = cellfun(@(x) x(1), NR_values);
R_values = cellfun(@(x) x(2), NR_values);

%% Generate multiscale filters

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

filters = cell(1, length(N_values)); %#ok<*NASGU>
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
    disp(['Working on filter at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    filters{i} = generateFourierWhitenFilter(images, 'NaturalImages', crtN, crtR);
end

%% Save multiscale filters

for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
    filter = filters{i};
    save(fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']), 'filter');
end

%% Load multiscale filters

filters = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
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

%% Generate the continuous stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_continuous = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on continuous stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_continuous{i} = analyzeImageSet(images, 'NaturalImages', crtN, ...
        'filter', filters{i}, 'equalize', 'contrast');
end

%% Save the continuous stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
res = res_continuous;
save(fullfile('save', 'natural_nosky_continuous_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Generate the binary stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_binary = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on binary stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_binary{i} = analyzeImageSet(images, 'NaturalImages', crtN, ...
        'filter', filters{i}, 'nLevels', 2, 'equalize', 'contrast');
end

%% Save the binary stats

filtersToString = cell2mat(arrayfun(@(n) ['_' int2str(n)], N_values, 'uniform', false));
res = res_binary;
save(fullfile('save', 'natural_nosky_binary_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Load binary and continuous stats

load(fullfile('save', 'natural_nosky_continuous_contrastadapt.mat'));
res_continuous = res;
cont_R_values = R_values;
cont_N_values = N_values;
cont_NR_values = NR_values;

load(fullfile('save', 'natural_nosky_binary_contrastadapt.mat'));
res_binary = res;
if ~isequal(R_values, cont_R_values) || ~isequal(N_values, cont_N_values) || ...
        ~isequal(NR_values, cont_NR_values)
    error('Inconsistent N, R sizes used between binary and continuous stats.');
end
clear('cont_R_values', 'cont_N_values', 'cont_NR_values', 'res');

%% Perform focus analysis

tic;
res_binary_with_focus = cell(size(res_binary));
for i = 1:length(res_binary_with_focus)
    disp(['Working on focus analysis for binary results ' int2str(i) ...
          '/' int2str(length(res_binary_with_focus)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_binary_with_focus{i} = runFocusAnalysis(res_binary{i});
end
disp(['Focus analysis for binary results took ' num2str(toc, '%.2f') ' seconds.']);

tic;
res_continuous_with_focus = cell(size(res_continuous));
for i = 1:length(res_continuous_with_focus)
    disp(['Working on focus analysis for continuous results ' int2str(i) ...
          '/' int2str(length(res_continuous_with_focus)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_continuous_with_focus{i} = runFocusAnalysis(res_continuous{i});
end
disp(['Focus analysis for continuous results took ' num2str(toc, '%.2f') ' seconds.']);

%% Save results with focus analysis

res = res_binary_with_focus;
save(fullfile('save', 'natural_nosky_binary_with_focus_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');

res = res_continuous_with_focus;
save(fullfile('save', 'natural_nosky_continuous_with_focus_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Load all stats

load(fullfile('save', 'natural_nosky_continuous_with_focus_contrastadapt.mat'));
res_continuous = res;
cont_R_values = R_values;
cont_N_values = N_values;
cont_NR_values = NR_values;

load(fullfile('save', 'natural_nosky_binary_with_focus_contrastadapt.mat'));
res_binary = res;
if ~isequal(R_values, cont_R_values) || ~isequal(N_values, cont_N_values) || ...
        ~isequal(NR_values, cont_NR_values)
    error('Inconsistent N, R sizes used between binary and continuous stats.');
end
clear('cont_R_values', 'cont_N_values', 'cont_NR_values', 'res');

% load continuous stats with histogram equalization
res_continuous_eq = getfield(open(fullfile('save', 'natural_nosky_continuous_with_focus.mat')), 'res');

%% Compare binary and continuous stats with contrast adaptation

for k = 1:length(NR_values)
    res0 = res_binary{k};
    res1 = res_continuous{k};
    
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
        
        smartscatter(res0.ev(:, i), res1.ev(:, i), 'alpha', 0.1, 'maxpoints', 50000);
        
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
    
    fig.Name = ['N=' int2str(N_values(k)) ', R=' int2str(R_values(k))];
    
    preparegraph;
end

%% Compare continuous stats with contrast adaptation vs. histogram equalization

for k = 1:length(NR_values)
    res0 = res_continuous{k};
    res1 = res_continuous_eq{k};
    
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
        
        smartscatter(res0.ev(:, i), res1.ev(:, i), 'alpha', 0.1, 'maxpoints', 50000);
        
        max_ax = max([abs(res0.ev(:, i)) ; abs(res1.ev(:, i))]);
        xlim([-max_ax max_ax]);
        ylim([-max_ax max_ax]);
        
        axis equal;
        xlabel('continuous, contrast adaptation');
        ylabel('continuous, histogram equalization');
        title(labels{i});
        beautifygraph;
        
        drawfitline(res0.ev(:, i), res1.ev(:, i), 'style', {'k', 'linewidth', 1});
    end
    
    fig.Name = ['N=' int2str(N_values(k)) ', R=' int2str(R_values(k))];
    
    preparegraph;
end

%% Generate the stochastic binary stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_stochastic_binary = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on stochastic binary stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_stochastic_binary{i} = analyzeImageSet(images, 'NaturalImages', crtN, ...
        'filter', filters{i}, 'nLevels', 2, 'quantizeType', 'stochastic', ...
        'equalize', 'contrast');
end

%% Perform focus analysis

tic;
res_stochastic_binary_with_focus = cell(size(res_stochastic_binary));
for i = 1:length(res_stochastic_binary_with_focus)
    disp(['Working on focus analysis for stochastic binary results ' int2str(i) ...
          '/' int2str(length(res_stochastic_binary_with_focus)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_stochastic_binary_with_focus{i} = runFocusAnalysis(res_stochastic_binary{i});
end
disp(['Focus analysis for stochastic binary results took ' num2str(toc, '%.2f') ' seconds.']);

%% Save the stochastic results with focus analysis

res = res_stochastic_binary_with_focus;
save(fullfile('save', 'natural_nosky_stochastic_binary_with_focus_contrastadapt.mat'), ...
    'res', 'R_values', 'N_values', 'NR_values');

%% Load stochastic binary stats

load(fullfile('save', 'natural_nosky_stochastic_binary_with_focus_contrastadapt.mat'));
res_stochastic_binary_with_focus = res;

clear('res');

%% Compare stochastic binary and continuous stats

%% ...make scatterplots for all 10 directions

NR_subset = 4;
% for k = 1:length(NR_values)
for k0 = 1:length(NR_subset)
    k = NR_subset(k0);
    
    res0 = res_stochastic_binary_with_focus{k};
    res1 = res_continuous{k};
    
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
        
        smartscatter(res0.ev(:, i), res1.ev(:, i), 'alpha', 0.1, 'maxpoints', 50000);
        
        max_ax = max([abs(res0.ev(:, i)) ; abs(res1.ev(:, i))]);
        xlim([-max_ax max_ax]);
        ylim([-max_ax max_ax]);
        
        axis equal;
        xlabel('stochastic binary');
        ylabel('continuous');
        title(labels{i});
        beautifygraph;
        
        drawfitline(res0.ev(:, i), res1.ev(:, i), 'style', {'k', 'linewidth', 1});
    end
    
    fig.Name = ['N=' int2str(N_values(k)) ', R=' int2str(R_values(k))];
    
    preparegraph;
end

%% Compare stochastic binary and binary stats

%% ...make scatterplots for all 10 directions

NR_subset = 4;
% for k = 1:length(NR_values)
for k0 = 1:length(NR_subset)
    k = NR_subset(k0);

    res0 = res_binary{k};
    res1 = res_stochastic_binary_with_focus{k};
    
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
        
        smartscatter(res0.ev(:, i), res1.ev(:, i), 'alpha', 0.1, 'maxpoints', 50000);
        
        max_ax = max([abs(res0.ev(:, i)) ; abs(res1.ev(:, i))]);
        xlim([-max_ax max_ax]);
        ylim([-max_ax max_ax]);
        
        axis equal;
        xlabel('binary');
        ylabel('stochastic binary');
        title(labels{i});
        beautifygraph;
        
        drawfitline(res0.ev(:, i), res1.ev(:, i), 'style', {'k', 'linewidth', 1});
    end
    
    fig.Name = ['N=' int2str(N_values(k)) ', R=' int2str(R_values(k))];
    
    preparegraph;
end