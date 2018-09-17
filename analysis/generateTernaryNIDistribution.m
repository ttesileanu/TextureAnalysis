% analyze natural images after ternarizing to generate a distribution over
% ternary textures

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   dbChoice
%       Choice of database to use. Available options are
%           'PennNoSky'     -- Penn image database, filtering out pictures
%                              with lots of sky
%           'vanHateren'    -- van Hateren image database, filtering out
%                              pictures with lots of sky or lots of
%                              human-made objects
%   valuesNR
%       Cell array of pairs [N, R] of downsampling factor and patch size
%       for which to run the analysis.
%   compressType
%       Choose the way in which image values are compressed in the [0, 1]
%       interval before ternarizing. Options are
%           'equalize' -- histogram equalization
%           'contrast' -- contrast adaptation

setdefault('dbChoice', 'PennNoSky');

setdefault('valuesNR', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});

setdefault('compressType', 'equalize');

%% Preprocess options

switch compressType
    case 'equalize'
        compressFunction = @(image) equalizeImage(image, 'jitter', 0, 'minLevels', 16);
    case 'contrast'
        compressFunction = @(image) contrastAdapt(image, 'minStd', 1e-3);
    otherwise
        error('Unknown compressType');
end

valuesN = cellfun(@(x) x(1), valuesNR);
valuesR = cellfun(@(x) x(2), valuesNR);

%% Load the filters

filters = cell(1, length(valuesN));
for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = valuesR(i);
    filterFilename = fullfile('filters', dbChoice, ...
        ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']);
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

%% Generate the ternary stats

images = parseImageNameFile([dbChoice 'Index.txt'], fullfile('NaturalImages', dbChoice));

results = cell(1, length(valuesN));
tic;
for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = size(filters{i}, 1);
    disp(['Working on ternary stats at N=' int2str(crtN) ', R=' int2str(crtR) ' (' int2str(i) ...
        '/' int2str(length(valuesN)) ')...']);
    
    % set up the preprocessing pipeline
    preprocessPipeline = {
        @logTransform, ...
        @(image) blockAverage(image, crtN), ...
        @(image) patchify(image, crtR), ...
        @(image) filterImage(image, filters{i}), ...
        compressFunction, ...
        @(image) quantizeImage(image, 3)};
    
    % image quantization has a small random component that matters for
    % patches that have many identical pixel values
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    results{i} = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipeline);
end
disp(['Generating the ternary stats took ' num2str(toc, '%.2f') ' seconds.']);

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

%% Save the ternary stats

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
save(fullfile('save', ['TernaryDistribution_' dbChoice compressExt '.mat']), 'results', ...
    'valuesR', 'valuesN', 'valuesNR');
