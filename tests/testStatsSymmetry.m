% some tests to make sure the symmetry calculations are working correctly

%% Load some data

% select preprocessing parameters
N = 2;
R = 32;

% load a filter
filter0 = open(fullfile('filters', sprintf('filter%dx%d.mat', N, R)));
filter = filter0.filter;

% select some images
images = parseImageNameFile('PennNoSkyIndex.txt', 'NaturalImages');
images = images(1:10);

%% Analyze textures without any transformation

% calculate texture distribution within selected images
preprocessPipeline = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3)};

% image quantization has a small random component that matters for
% patches that have many identical pixel values
% to keep things reproducible, we fix the random number generator seed
rng('default');
results = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipeline);

%% Analyze textures after LR flip

preprocessPipelineLR = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) deal(fliplr(patches), []), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3)};

rng('default');
resultsLR = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineLR);

%% Compare

[fastEvLR, shuffleLR] = applyStatsGeometricPermutation(results.ev, 3, 'BADC');
scatterfit(fastEvLR, resultsLR.ev);

%% Analyze textures after UD flip

preprocessPipelineUD = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) deal(flipud(patches), []), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3)};

rng('default');
resultsUD = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineUD);

%% Compare

[fastEvUD, shuffleUD] = applyStatsGeometricPermutation(results.ev, 3, 'CDAB');
scatterfit(fastEvUD, resultsUD.ev);

%% Analyze textures after 90 degree CW rotation

preprocessPipelineRotCW = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) deal(rot90(patches, -1), []), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3)};

rng('default');
resultsRotCW = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineRotCW);

%% Compare

[fastEvRotCW, shuffleRotCW] = applyStatsGeometricPermutation(results.ev, 3, 'BDAC');
scatterfit(fastEvRotCW, resultsRotCW.ev);

%% Analyze textures after color cycle

preprocessPipelineColorCycle = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3), ...
    @(patches) deal(applyImageColorTransformation(patches, 3, 1, 1), [])};

rng('default');
resultsColorCycle = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineColorCycle);

%% Compare

[fastEvColorCycle, shuffleColorCycle] = applyStatsColorTransformation(results.ev, 3, 1, 1);
scatterfit(fastEvColorCycle, resultsColorCycle.ev);

%% Analyze textures after gray-white flip

preprocessPipelineGWFlip = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3), ...
    @(patches) deal(applyImageColorTransformation(patches, 3, 2, 0), [])};

rng('default');
resultsGWFlip = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineGWFlip);

%% Compare

[fastEvGWFlip, shuffleGWFlip] = applyStatsColorTransformation(results.ev, 3, 2, 0);
scatterfit(fastEvGWFlip, resultsGWFlip.ev);

%% Analyze textures after black-gray flip

preprocessPipelineBGFlip = {
    @logTransform, ...
    @(image) blockAverage(image, N), ...
    @(image) patchify(image, R), ...
    @(patches) filterImage(patches, filter), ...
    @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
    @(patches) quantizeImage(patches, 3), ...
    @(patches) deal(applyImageColorTransformation(patches, 3, 2, 1), [])};

rng('default');
resultsBGFlip = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineBGFlip);

%% Compare

[fastEvBGFlip, shuffleBGFlip] = applyStatsColorTransformation(results.ev, 3, 2, 1);
scatterfit(fastEvBGFlip, resultsBGFlip.ev);

%% Analyze textures after all color permutations

resultsColorPermutation = cell(2, 3);
for xCoeff = 1:2
    for yCoeff = 0:2
        preprocessPipelineColorPermutation = {
            @logTransform, ...
            @(image) blockAverage(image, N), ...
            @(image) patchify(image, R), ...
            @(patches) filterImage(patches, filter), ...
            @(patches) equalizeImage(patches, 'jitter', 0, 'minLevels', 16), ...
            @(patches) quantizeImage(patches, 3), ...
            @(patches) deal(applyImageColorTransformation(patches, 3, xCoeff, yCoeff), [])};
        rng('default');
        resultsColorPermutation{xCoeff, yCoeff+1} = ...
            generateTextureDistribution(images, 3, 'preprocessing', preprocessPipelineColorPermutation);
    end
end

%% Compare

for xCoeff = 1:2
    for yCoeff = 0:2        
        [fastEvColorPermutation, shuffleColorPermutation] = ...
            applyStatsColorTransformation(results.ev, 3, xCoeff, yCoeff);
        subplot(2, 3, (xCoeff-1)*3 + yCoeff + 1);
        scatterfit(fastEvColorPermutation, resultsColorPermutation{xCoeff, yCoeff+1}.ev);
    end
end