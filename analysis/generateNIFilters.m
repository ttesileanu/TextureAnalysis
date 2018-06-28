% generate the whitening filters from natural images

%% Select downsampling factor (N) and patch size (R)

% [N, R] pairs
valuesNR = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]};
% valuesNR = {[1, 32]};

valuesN = cellfun(@(x) x(1), valuesNR);
valuesR = cellfun(@(x) x(2), valuesNR);

%% Generate filters for every choice of parameters

images = parseImageNameFile('PennNoSkyIndex.txt', 'NaturalImages');

filters = cell(1, length(valuesN)); %#ok<*NASGU>
tic;
for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = valuesR(i);
    disp(['Working on filter at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(valuesN)) ')...']);
    
    % set up the preprocessing pipeline
    preprocessPipeline = {
        @logTransform, ...
        @(image) blockAverage(image, crtN), ...
        @(image) patchify(image, crtR)};
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    filters{i} = generateFourierWhitenFilter(images, 'preprocessing', preprocessPipeline);
end
disp(['Filter generation took ' num2str(toc, '%.2f') ' seconds.']);

%% Save the filters

for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = valuesR(i);
    filter = filters{i};
    save(fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']), 'filter');
end