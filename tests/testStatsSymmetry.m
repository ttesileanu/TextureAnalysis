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