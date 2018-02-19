% generate results from G=3 analysis

%% Select blockAF and patch sizes

% [N, R] pairs
% NR_values = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
%     [4, 32], [4, 48], [4, 64]};
NR_values = {[2, 32]};

N_values = cellfun(@(x) x(1), NR_values);
R_values = cellfun(@(x) x(2), NR_values);

%% Generate multiscale non-log filters

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
    filters{i} = generateFourierWhitenFilter(images, 'NaturalImages', crtN, crtR, ...
        'doLog', false);
end

%% Save multiscale non-log filters

for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
    filter = filters{i};
    save(fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '_nolog.mat']), 'filter');
end

%% Load multiscale filters

filters = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = R_values(i);
    filterFilename = fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '_nolog.mat']);
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

%% Generate the ternary stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');

res_ternary = cell(1, length(N_values));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    disp(['Working on ternary stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
%     res_ternary{i} = analyzeImageSet(images, 'NaturalImages', crtN, 'nLevels', 3, ...
%         'filter', filters{i}, 'doLog', false);
    res_ternary{i} = analyzeImageSet(images, 'NaturalImages', crtN, 'nLevels', 3, ...
        'filter', filters{i}, 'doLog', false, 'equalize', 'contrast');
end

%% Save the ternary stats

res = res_ternary;
% save(fullfile('save', 'natural_nosky_ternary_nolog.mat'), 'res', 'R_values', 'N_values', 'NR_values');
save(fullfile('save', 'natural_nosky_ternary_nolog_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');

%% Load ternary stats

% load(fullfile('save', 'natural_nosky_ternary_nolog.mat'));
load(fullfile('save', 'natural_nosky_ternary_nolog_contrastadapt.mat'));
res_ternary = res;
clear('res');

%% Compare to John's results

john = open(fullfile('/Users/ttesileanu/Dropbox/Textures/GrayLevelNIStats/RawNIG3Stats', ...
    'NaturalImageStatsNSL_G3_N2_PC1_32_32.mat'));

% reduce John's data to only the two components that we keep
john.stats_red = reshape(permute(john.stats(:, 1:2, :), [3 2 1]), [], size(res_ternary{1}.ev, 2));

% need to reorganize patches because John's data goes in a different order
res_reorg = transposeImagePatches(res_ternary{1});

% differences are due to the following
% * John's code uses periodic b.c. for estimating statistics; this code doesn't
% * John's equalization and quantization routine gives slightly different
%   results than the one used here
% * there probably are differences for patches that have a lot of identical
%   pixels, which in this code get randomized

scatterfit(john.stats_red, res_reorg.ev, 'fitopts', {'line', [1, 0]});
xlabel('John analysis');
ylabel('This analysis');
beautifygraph;
preparegraph;