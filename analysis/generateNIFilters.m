% generate the whitening filters from natural images

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

setdefault('dbChoice', 'PennNoSky');

setdefault('valuesNR', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});

%% Preprocess options

valuesN = cellfun(@(x) x(1), valuesNR);
valuesR = cellfun(@(x) x(2), valuesNR);

%% Generate filters for every choice of parameters

images = parseImageNameFile([dbChoice 'Index.txt'], fullfile('NaturalImages', dbChoice));

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
    
    % XXX random number shouldn't be used here, so this is overkill
    rng('default');
    filters{i} = generateFourierWhitenFilter(images, 'preprocessing', preprocessPipeline);
end
disp(['Filter generation took ' num2str(toc, '%.2f') ' seconds.']);

%% Save the filters

for i = 1:length(valuesN)
    crtN = valuesN(i);
    crtR = valuesR(i);
    filter = filters{i};
    save(fullfile('filters', dbChoice, ['filter' int2str(crtN) 'x' int2str(crtR) '.mat']), 'filter');
end