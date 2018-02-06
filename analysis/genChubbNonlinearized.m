% generate continuous stats for images that have been passed through
% Chubb-like nonlinearities

%% Select blockAF and patch sizes

% [N, R] pairs
% NR_values = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
%     [4, 32], [4, 48], [4, 64]};
NR_values = {[2, 32]};

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

%% Generate the nonlinearized continuous stats

images = parseImageNameFile('Natural_Images_Large_NoSky_Index.txt', 'NaturalImages');
% images = parseImageNameFile('Natural_Images_Test_Index.txt', 'NaturalImages');

chubb_nonlin = open('data/chubb_nonlinearities.mat');

% first and second nonlinearities are identical
chubb_mask = [2 3 4];
chubb_nonlin.pvals = chubb_nonlin.pvals(:, chubb_mask);

res_nonlinear = cell(length(N_values), size(chubb_nonlin.pvals, 2));
for i = 1:length(N_values)
    crtN = N_values(i);
    crtR = size(filters{i}, 1);
    
    for j = 1:size(chubb_nonlin.pvals, 2)
        disp(['Working on nonlinearized stats at N=' int2str(crtN) ', R=' int2str(crtR) ...
            ', nonlin=' int2str(j) ' (' int2str(i) ...
        '/' int2str(length(N_values)) ')...']);
        
        % image quantization has a small random component that matters for
        % patches that have many identical pixel values
        % to keep things reproducible, we fix the random number generator seed
        rng('default');
        res_nonlinear{i, j} = analyzeImageSet(images, 'NaturalImages', crtN, ...
            'filter', filters{i}, 'nonlinearity', chubb_nonlin.pvals(:, j), ...
            'equalize', 'contrast');
        res_nonlinear{i, j}.nonlinearity = chubb_nonlin.pvals(:, j);
    end
end

%% Save the nonlinearized continuous stats

res = res_nonlinear;
% save(fullfile('save', 'natural_nosky_nonlinear.mat'), 'res', 'R_values', 'N_values', 'NR_values');
save(fullfile('save', 'natural_nosky_nonlinear_contrastadapt.mat'), 'res', 'R_values', 'N_values', 'NR_values');
clear('res');