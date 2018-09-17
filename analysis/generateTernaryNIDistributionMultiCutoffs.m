% analyze natural images after ternarizing to generate a distribution over
% ternary textures

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   valuesNR
%       Cell array of pairs [N, R] of downsampling factor and patch size
%       for which to run the analysis.
%   cutoffChoices
%       Choices to use for the quantization cutoffs.

setdefault('valuesNR', {[2, 32]});

defaultCutoffs0 = arrayfun(@(g) [(1-g)/2, (1+g)/2], linspace(0, 1, 25), 'uniform', false);
defaultCutoffs = cutoffs0(2:end-1);
setdefault('cutoffChoices', defaultCutoffs);

%% Preprocess options

valuesN = cellfun(@(x) x(1), valuesNR);
valuesR = cellfun(@(x) x(2), valuesNR);

%% Load the filters

filters = cell(1, length(valuesN));
for i = 1:length(valuesN)
    N = valuesN(i);
    R = valuesR(i);
    filterFilename = fullfile('filters', ['filter' int2str(N) 'x' int2str(R) '.mat']);
    % load the filter from whichever variable has the right size and type
    crtFilter = open(filterFilename);
    fields = fieldnames(crtFilter);
    isfilter = @(f) isnumeric(f) && ismatrix(f) && all(size(f) == [R R]);
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

images = parseImageNameFile('PennNoSkyIndex.txt', fullfile('NaturalImages', 'PennNoSky'));

tic;
for i = 1:length(valuesN)
    N = valuesN(i);
    R = size(filters{i}, 1);
    disp(['Working on ternary stats at N=' int2str(N) ', R=' int2str(R) ' (' int2str(i) ...
        '/' int2str(length(valuesN)) ')...']);
    
    sharpness = [];
    for k = 1:length(cutoffChoices)
        crtCutoffs = cutoffChoices{k};
        
        disp(['...cutoffs [' num2str(crtCutoffs, '%6.2f') '], time = ' num2str(toc, '%.2f') ' seconds...']);
        
        % set up the preprocessing pipeline
        preprocessPipeline = {
            @logTransform, ...
            @(image) blockAverage(image, N), ...
            @(image) patchify(image, R), ...
            @(image) filterImage(image, filters{i}), ...
            @(image) equalizeImage(image, 'jitter', 0, 'minLevels', 64), ...
            @(image) quantizeImage(image, 3, 'cutoffs', crtCutoffs)};
        
        % image quantization has a small random component that matters for
        % patches that have many identical pixel values
        % to keep things reproducible, we fix the random number generator seed
        rng('default');
        results = generateTextureDistribution(images, 3, 'preprocessing', preprocessPipeline);
        
        % perform focus analysis
        disp(['Working on focus analysis for ternary results ' int2str(i) ...
            '/' int2str(length(results)) ', cutoffs [' num2str(crtCutoffs, '%6.2f') ']...']);
        % focus analysis has a random component
        % to keep things reproducible, we fix the random number generator seed
        rng('default');
        try
            if isempty(sharpness)
                results = runFocusAnalysis(results);
                sharpness = results.focus.sharpness;
            else
                results.focus.sharpness = sharpness;
                results = runFocusAnalysis(results, 'reuseSharpness', true);
            end
        catch me
            warning(['Focus analysis failed, message: ' me.message]);
        end
        
        % save
        disp('Saving to file...');
        filename = ['TernaryDistribution_PennNoSky_' int2str(N) 'x' int2str(R) ...
            '_cutoff_' int2str(round(100*crtCutoffs(1))) '_' int2str(round(100*crtCutoffs(2))) '.mat'];
        
        cutoffs = crtCutoffs;
        
        save(fullfile('save', 'multicutoff', filename), 'results', 'N', 'R', 'cutoffs');
    end
end
disp(['Generating the ternary stats took ' num2str(toc, '%.2f') ' seconds.']);