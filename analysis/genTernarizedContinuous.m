% generate natural statistics of textures using continuous analysis, but
% after ternarizing the natural patches

%% Select blockAF and patch sizes

% [N, R] pairs
NR_values = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]};

N_values = cellfun(@(x) x(1), NR_values);
R_values = cellfun(@(x) x(2), NR_values);

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

%% Generate the continuous stats of ternarized images

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_ternarized_continuous = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on continuous stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_ternarized_continuous{i} = analyzeImageSet(images, 'NaturalImages', crtN, ...
        'filter', filters{i}, 'quantize', 3, 'nLevels', inf);
end

%% Perform focus analysis

tic;
res_ternarized_continuous_with_focus = cell(size(res_ternarized_continuous));
for i = 1:length(res_ternarized_continuous_with_focus)
    disp(['Working on focus analysis for continuous results ' int2str(i) ...
          '/' int2str(length(res_ternarized_continuous_with_focus)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    res_ternarized_continuous_with_focus{i} = runFocusAnalysis(res_ternarized_continuous{i});
end
disp(['Focus analysis took ' num2str(toc, '%.2f') ' seconds.']);

%% Save results with focus analysis

res = res_ternarized_continuous_with_focus;
save(fullfile('save', 'natural_nosky_ternarized_continuous_with_focus.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Load all stats

load(fullfile('save', 'natural_nosky_ternarized_continuous_with_focus.mat'));
res_ternarized_continuous = res;
cont_R_values = R_values;
cont_N_values = N_values;
cont_NR_values = NR_values;

load(fullfile('save', 'natural_nosky_continuous_with_focus.mat'));
res_continuous = res;
if ~isequal(R_values, cont_R_values) || ~isequal(N_values, cont_N_values) || ...
        ~isequal(NR_values, cont_NR_values)
    error('Inconsistent N, R sizes used between ternarized and non-ternarized continuous stats.');
end
clear('cont_R_values', 'cont_N_values', 'cont_NR_values', 'res');

%% Compare ternarized and non-ternarized stats

%% ...make scatterplots for all 10 directions

for k = 1:length(NR_values)
    res0 = res_continuous{k};
    res1 = res_ternarized_continuous{k};
    
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
        xlabel('non-ternarized');
        ylabel('ternarized');
        title(labels{i});
        beautifygraph;
        
        drawfitline(res0.ev(:, i), res1.ev(:, i), 'style', {'k', 'linewidth', 1});
    end
    
    fig.Name = ['N=' int2str(N_values(k)) ', R=' int2str(R_values(k))];
    
    preparegraph;
end