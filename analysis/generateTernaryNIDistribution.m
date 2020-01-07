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
%   filterScope
%       Choose the scope of the whitening filter:
%           'patch' (default) -- whiten each patch separately
%           'image'           -- whiten before patchifying
%   compressScope
%       Choose the scope of the compression. Note that setting this to
%       'image' forces `filterScope` to be 'image' as well.
%           'patch' (default) -- compress each patch separately
%           'image'           -- compress before patchifying

setdefault('dbChoice', 'PennNoSky');

setdefault('valuesNR', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});

setdefault('compressType', 'equalize');

setdefault('compressScope', 'patch');

setdefault('filterScope', 'patch');

%% Preprocess options

switch compressType
    case 'equalize'
        compressFunction = @(image) equalizeImage(image, 'jitter', 0, 'minLevels', 16);
    case 'contrast'
        compressFunction = @(image) contrastAdapt(image, 'minStd', 1e-3);
    otherwise
        error('Unknown compressType');
end

if strcmp(compressScope, 'image')
    filterScope = 'image';
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
        @(image) blockAverage(image, crtN)};
    
    if strcmp(filterScope, 'image')
        preprocessPipeline = [preprocessPipeline {@(image) filterImage(image, filters{i})}]; %#ok<AGROW>
    end
    if strcmp(compressScope, 'image')
        preprocessPipeline = [preprocessPipeline {compressFunction}]; %#ok<AGROW>
    end
    
    preprocessPipeline = [preprocessPipeline {@(image) patchify(image, crtR)}]; %#ok<AGROW>
    
    if strcmp(filterScope, 'patch')
        preprocessPipeline = [preprocessPipeline {@(image) filterImage(image, filters{i})}]; %#ok<AGROW>
    end
    if strcmp(compressScope, 'patch')
        preprocessPipeline = [preprocessPipeline {compressFunction}]; %#ok<AGROW>
    end
    
    preprocessPipeline = [preprocessPipeline {@(image) quantizeImage(image, 3)}]; %#ok<AGROW>
    
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
    disp(['Working on focus analysis for ternary results ' int2str(i) ...
          '/' int2str(length(results)) '...']);
    % focus analysis has a random component
    % to keep things reproducible, we fix the random number generator seed
    rng('default');
    resultsFocus{i} = runFocusAnalysis(results{i});
end
disp(['Focus analysis for ternary stats took ' num2str(toc, '%.2f') ' seconds.']);
results = resultsFocus;

%% Save the ternary stats

if ~strcmp(compressType, 'equalize')
    extras = ['_' compressType];
else
    extras = '';
end
if ~strcmp(filterScope, 'patch')
    extras = [extras '_flt' filterScope];
end
if ~strcmp(compressScope, 'patch')
    extras = [extras '_comp' compressScope];
end
filename = ['TernaryDistribution_' dbChoice extras '.mat'];
save(fullfile('save', filename), 'results', ...
    'valuesR', 'valuesN', 'valuesNR', '-v7.3');
