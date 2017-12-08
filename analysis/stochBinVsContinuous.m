% compare averaged results from binary statistics applied to stochastically
% binarized image patches, to results from the continuous stats

%% Load an image, to have a good source of patches

% make this reproducible -- patches that are very uniform might lead to
% stochastic results from the equalization
rng(1);

filter = open(fullfile('filters', 'filter1x32.mat'));
image = preprocessImage(fullfile('NaturalImages', 'cd01A', 'DSC_0002_LUM.mat'), 1, ...
    'filter', filter.filter, 'equalize', true);

%% Run continuous stats analysis on patches from this image

patchSize = 32;
continuous_res = analyzePatches(image, inf, patchSize);

%% Run many samples of binary stats on stochastic binarization

% make this reproducible
rng(3487);

nSamples = 250;
stochastic_binary_res = cell(1, nSamples);

tic;
n_max_dots = 40;
last_dots = 0;
fprintf('Progress: ');
for i = 1:nSamples
    crt_dots = round((i-1)/(nSamples-1)*n_max_dots);
    if crt_dots > last_dots
        fprintf(1, repmat('.', 1, crt_dots - last_dots));
        last_dots = crt_dots;
    end
    crt_image = stochasticBinarize(image);
    
    stochastic_binary_res{i} = analyzePatches(crt_image, 2, patchSize);
end
fprintf('\n');
disp(['Calculation took ' num2str(toc, '%.2f') ' seconds.']);

%% Calculate averages of stochastic binary stats

cell_stochastic_binary_ev = cellfun(@(s) s.ev, stochastic_binary_res, 'uniform', false);
array_stochastic_binary_ev = cat(3, cell_stochastic_binary_ev{:});
mean_stochastic_binary_ev = mean(array_stochastic_binary_ev, 3);

%% Compare

scatterfit(continuous_res.ev, mean_stochastic_binary_ev);