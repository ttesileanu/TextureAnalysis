% analyze natural images after binarizing to generate a distribution over
% binary textures

%% Select downsampling factor (N) and patch size (R)

% [N, R] pairs
% valuesNR = {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
%     [4, 32], [4, 48], [4, 64]};
valuesNR = {[1, 32]};

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

%% Load the filters

filters = cell(1, length(valuesN));
for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = valuesR(i);
    filterFilename = fullfile('filters', ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']);
    % load the filter from whichever variable has the right size and type
    crtFilter = open(filterFilename);
    fields = fieldnames(crtFilter);
    isfilter = @(f) isnumeric(f) && ismatrix(f) && all(size(f) == [crtR crtR]);
    valid = cellfun(@(s) isfilter(crtFilter.(s)), fields);
    if sum(valid) == 0
        error(['Can''t find valid filter data in ' filterFilename '.']);
    elseif sum(valid) > 1
        error(['Don''t know which variable to use from ' filterFilename '.']);
    end
    filterField = fields{valid};
    filters{i} = crtFilter.(filterField);
end

%% Generate the binary stats

images = parseImageNameFile('PennNoSkyIndex.txt', 'NaturalImages');

results = cell(1, length(valuesN));
tic;
for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = size(filters{i}, 1);
    disp(['Working on binary stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(valuesN)) ')...']);
    
    % set up the preprocessing pipeline
    preprocessPipeline = {
        @logTransform, ...
        @(image) blockAverage(image, crtN), ...
        @(image) patchify(image, crtR), ...
        @(image) filterImage(image, filters{i}), ...
        @equalizeImage, ...
        @(image) quantizeImage(image, 2)};
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    results{i} = generateTextureDistribution(images, 2, 'preprocessing', preprocessPipeline);
end
disp(['Generating the binary stats took ' num2str(toc, '%.2f') ' seconds.']);

%% Perform focus analysis

tic;
resultsFocus = cell(1, length(valuesN));
for i = 1:length(results)
    disp(['Working on focus analysis for binary results ' int2str(i) ...
          '/' int2str(length(results)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    resultsFocus{i} = runFocusAnalysis(results{i});
end
disp(['Focus analysis for binary stats took ' num2str(toc, '%.2f') ' seconds.']);
results = resultsFocus;

%% Save the binary stats

save(fullfile('save', 'BinaryDistribution_PennNoSky.mat'), 'results', 'R_values', 'N_values', 'NR_values');