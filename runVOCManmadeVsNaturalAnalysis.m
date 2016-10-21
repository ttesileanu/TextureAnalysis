function stats = runVOCManmadeVsNaturalAnalysis(res, varargin)
% runVOCManmadeVsNaturalAnalysis Run an analysis comparing texture
% statistics for man-made vs. natural objects in the VOC database.
%
%   Options:
%    'kfold': integer
%       Number of folds to use for k-fold cross-validation.
%    'nreps': integer
%       Number of times to repeat the cross-validation (each time with a
%       different random subsampling).
%    'skipEv1': logical
%       If true (the default), skip the first component of the ev vector
%       (the luminance, gamma) in the analyses. Set to false to use all the
%       ev data.
%    'returnObj': logical
%       If true, return the Gaussian mixture objects returned by the linear
%       discrimination analysis.
%    'progressEvery': double
%       How often to display progress information (in seconds), after the
%       'progressStart' period (see below) elapsed.
%    'progressStart': double
%       How long to wait before displaying progress information for the
%       first time. Set to infinity to never display progress.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('kfold', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('nreps', 100, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('skipEv1', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('returnObj', false, @(b) islogical(b) && isscalar(b));

parser.addParameter('progressEvery', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('progressStart', 20, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% find the mask for all man-made and natural patches
allManmade = ismember(res.objIds, [2 4]);
allNatural = ismember(res.objIds, 1);
allManmadeIdxs = find(allManmade);
allNaturalIdxs = find(allNatural);
nManmade = length(allManmadeIdxs);
nNatural = length(allNaturalIdxs);
disp(['The dataset contains ' int2str(nManmade) ' man-made patches and ' ...
    int2str(nNatural) ' natural patches.']);

% number of patches to analyze is the minimum over patch types
% adjusted so that it's divisible by the number of folds
nPatches = floor(min(nManmade, nNatural) / params.kfold)*params.kfold;
disp(['Using ' int2str(nPatches) ' patches of each type.']);

% repeat the analysis params.nreps times, using a different subsampling
% each time
ev = res.ev(:, (1+params.skipEv1):end);
stats.raw.natural.allPC = zeros(size(ev, 2), size(ev, 2), params.nreps, params.kfold);
stats.raw.manmade.allPC = zeros(size(ev, 2), size(ev, 2), params.nreps, params.kfold);
stats.raw.natural.allSVals = zeros(size(ev, 2), params.nreps, params.kfold);
stats.raw.manmade.allSVals = zeros(size(ev, 2), params.nreps, params.kfold);
stats.raw.lda.obj = cell(params.nreps, params.kfold);
stats.raw.lda.trainingPerf = zeros(params.nreps, params.kfold);
stats.raw.lda.validationPerf = zeros(params.nreps, params.kfold);
stats.raw.lda.trainingContingency = zeros(2, 2, params.nreps, params.kfold);
stats.raw.lda.validationContingency = zeros(2, 2, params.nreps, params.kfold);
stats.raw.lda.naturalMu = zeros(size(ev, 2), params.nreps, params.kfold);
stats.raw.lda.manmadeMu = zeros(size(ev, 2), params.nreps, params.kfold);
stats.raw.lda.discriminant = zeros(size(ev, 2), params.nreps, params.kfold);
stats.raw.lda.const = zeros(params.nreps, params.kfold);

t0 = tic;
tEvery = tic;
progressWritten = false;

for i = 1:params.nreps
    % output progress information if required
    if (~progressWritten && toc(t0) > params.progressStart) || ...
            (progressWritten && toc(tEvery) > params.progressEvery)
        disp(['Working on repetition ' int2str(i) ' of ' int2str(params.nreps) ...
            ', elapsed ' num2str(toc(t0), '%.1f') ' seconds...']);
        progressWritten = true;
        tEvery = tic;
    end
    
    % figure out the folds
    naturalOrdering = allNaturalIdxs(randperm(nNatural));
    manmadeOrdering = allManmadeIdxs(randperm(nManmade));
    
    naturalFolds = reshape(naturalOrdering(1:nPatches), [], params.kfold);
    manmadeFolds = reshape(manmadeOrdering(1:nPatches), [], params.kfold);
    
    % do the analysis keeping each fold as a validation set
    for k = 1:params.kfold
        % figure out the training and validation datasets
        naturalTrainingIdxs = flatten(naturalFolds(:, setdiff(1:params.kfold, k)));
        manmadeTrainingIdxs = flatten(manmadeFolds(:, setdiff(1:params.kfold, k)));
        
        naturalValidationIdxs = naturalFolds(:, k);
        manmadeValidationIdxs = manmadeFolds(:, k);
        
        naturalTrainingEv = ev(naturalTrainingIdxs, :);
        manmadeTrainingEv = ev(manmadeTrainingIdxs, :);
        
        naturalValidationEv = ev(naturalValidationIdxs, :);
        manmadeValidationEv = ev(manmadeValidationIdxs, :);
        
        % get principal components
        [stats.raw.natural.allPC(:, :, i, k), stats.raw.natural.allSVals(:, i, k)] = ...
            doPca(naturalTrainingEv);
        [stats.raw.manmade.allPC(:, :, i, k), stats.raw.manmade.allSVals(:, i, k)] = ...
            doPca(manmadeTrainingEv);
        
        trainingEv = [naturalTrainingEv ; manmadeTrainingEv];
        validationEv = [naturalValidationEv ; manmadeValidationEv];
        [stats.raw.both.allPC(:, :, i, k), stats.raw.both.allSVals(:, i, k)] = ...
            doPca(trainingEv);
        
        % perform linear discriminant analysis
        % label 0 = natural; label 1 = man-made
        trainingLabels = [zeros(size(naturalTrainingEv, 1), 1) ; ones(size(manmadeTrainingEv, 1), 1)];
        validationLabels = [zeros(size(naturalValidationEv, 1), 1) ; ones(size(manmadeValidationEv, 1), 1)];
        stats.raw.lda.obj{i, k} = fitcdiscr(trainingEv, trainingLabels);
        
        % check performance of LDA on both training and validation datasets
        ldaPredsTraining = stats.raw.lda.obj{i, k}.predict(trainingEv);
        ldaPredsValidation = stats.raw.lda.obj{i, k}.predict(validationEv);
        stats.raw.lda.trainingPerf(i, k) = mean(ldaPredsTraining == trainingLabels);
        stats.raw.lda.validationPerf(i, k) = mean(ldaPredsValidation == validationLabels);
        stats.raw.lda.trainingContingency(:, :, i, k) = getContingencyTable(ldaPredsTraining, trainingLabels);
        stats.raw.lda.validationContingency(:, :, i, k) = getContingencyTable(ldaPredsValidation, validationLabels);
        naturalLdaCoeff = find(stats.raw.lda.obj{i, k}.ClassNames == 0);
        manmadeLdaCoeff = find(stats.raw.lda.obj{i, k}.ClassNames == 1);
        stats.raw.lda.naturalMu(:, i, k) = stats.raw.lda.obj{i, k}.Mu(naturalLdaCoeff);
        stats.raw.lda.manmadeMu(:, i, k) = stats.raw.lda.obj{i, k}.Mu(manmadeLdaCoeff);
        stats.raw.lda.discriminant(:, i, k) = ...
            stats.raw.lda.obj{i, k}.Coeffs(naturalLdaCoeff, manmadeLdaCoeff).Linear;
        stats.raw.lda.const(i, k) = ...
            stats.raw.lda.obj{i, k}.Coeffs(naturalLdaCoeff, manmadeLdaCoeff).Const;
    end
end

if progressWritten
    disp(['Analysis finished, took ' num2str(toc(t0), '%.2f') ' seconds.']);
end

% summarize principal component data
fields = {'natural', 'manmade', 'both'};
for i = 1:length(fields)
    crtField = fields{i};
    
    % summarize data per repetition by averaging over k-folds
    stats.reps.(crtField).allPC = mean(stats.raw.(crtField).allPC, 4);
    stats.reps.(crtField).allSVals = mean(stats.raw.(crtField).allSVals, 3);
    
    % summarize all the analysis results
    stats.summary.(crtField).allPC = mean(stats.reps.(crtField).allPC, 3);
    stats.summary.(crtField).allSVals = mean(stats.reps.(crtField).allSVals, 2);
    
    % include standard deviations in the summaries
    stats.summary.(crtField).stdAllPC = std(stats.reps.(crtField).allPC, [], 3);
    stats.summary.(crtField).stdAllSVals = std(stats.reps.(crtField).allSVals, [], 2);
end

% get top principal component, for convenience
for i = 1:length(fields)
    crtField = fields{i};
    
    stats.raw.(crtField).topPC = squeeze(stats.raw.(crtField).allPC(:, 1, :, :));
    stats.reps.(crtField).topPC = squeeze(stats.reps.(crtField).allPC(:, 1, :));
    stats.summary.(crtField).topPC = stats.summary.(crtField).allPC(:, 1);
    stats.summary.(crtField).stdTopPC = stats.summary.(crtField).stdAllPC(:, 1);
end

% summarize data per repetition by averaging over k-folds
stats.reps.lda.trainingPerf = mean(stats.raw.lda.trainingPerf, 2);
stats.reps.lda.validationPerf = mean(stats.raw.lda.validationPerf, 2);
stats.reps.lda.trainingContingency = mean(stats.raw.lda.trainingContingency, 4);
stats.reps.lda.validationContingency = mean(stats.raw.lda.validationContingency, 4);
stats.reps.lda.naturalMu = mean(stats.raw.lda.naturalMu, 3);
stats.reps.lda.manmadeMu = mean(stats.raw.lda.manmadeMu, 3);
stats.reps.lda.discriminant = mean(stats.raw.lda.discriminant, 3);
stats.reps.lda.const = mean(stats.raw.lda.const, 2);

% summarize all the analysis results
stats.summary.lda.trainingPerf = mean(stats.reps.lda.trainingPerf);
stats.summary.lda.validationPerf = mean(stats.reps.lda.validationPerf);
stats.summary.lda.trainingContingency = mean(stats.reps.lda.trainingContingency, 3);
stats.summary.lda.validationContingency = mean(stats.reps.lda.validationContingency, 3);
stats.summary.lda.naturalMu = mean(stats.reps.lda.naturalMu, 2);
stats.summary.lda.manmadeMu = mean(stats.reps.lda.manmadeMu, 2);
stats.summary.lda.discriminant = mean(stats.reps.lda.discriminant, 2);
stats.summary.lda.const = mean(stats.reps.lda.const);

% include standard deviations in the summaries
stats.summary.lda.stdTrainingPerf = std(stats.reps.lda.trainingPerf);
stats.summary.lda.stdValidationPerf = std(stats.reps.lda.validationPerf);
stats.summary.lda.stdTrainingContingency = std(stats.reps.lda.trainingContingency, [], 3);
stats.summary.lda.stdValidationContingency = std(stats.reps.lda.validationContingency, [], 3);
stats.summary.lda.stdNaturalMu = std(stats.reps.lda.naturalMu, [], 2);
stats.summary.lda.stdManmadeMu = std(stats.reps.lda.manmadeMu, [], 2);
stats.summary.lda.stdDiscriminant = std(stats.reps.lda.discriminant, [], 2);
stats.summary.lda.stdConst = std(stats.reps.lda.const);

% cosine similarities between various vectors
% order: PC1 both datasets, PC1 natural, PC1 manmade, discriminant
stats.reps.similarity.matrix = zeros(4, 4, params.nreps);
for i = 1:params.nreps
    crtVecs = [stats.reps.both.topPC(:, i) ...
               stats.reps.natural.topPC(:, i) ...
               stats.reps.manmade.topPC(:, i) ...
               stats.reps.lda.discriminant(:, i) / norm(stats.reps.lda.discriminant(:, i))];
    stats.reps.similarity.matrix(:, :, i) = crtVecs'*crtVecs;
end
stats.summary.similarity.matrix = mean(stats.reps.similarity.matrix, 3);
stats.summary.similarity.stdMatrix = std(stats.reps.similarity.matrix, [], 3);

if ~params.returnObj
    % get rid of LDA objects
    % XXX would be better to not save these in the first place
    stats.raw.lda = rmfield(stats.raw.lda, 'obj');
end

end

function [V, S, mu] = doPca(ev)
% Return the right singular vectors and the singular values for the data.
%
% The signs of the singular vectors are chosen so that the entry that is
% largest in absolute value has positive sign.

mu = mean(ev, 1);
evCentered = bsxfun(@minus, ev, mu);

[~, S, V] = svd(evCentered, 0);

S = diag(S);

for i = 1:size(V, 2)
    [~, idx] = max(abs(V(:, i)));
    if V(idx, i) < 0
        V(:, i) = -V(:, i);
    end
end

end

function M = getContingencyTable(preds, truths)
% Return a 2x2 contingency table showing how often predictions match the
% ground truth. "Truth" indices are along the rows, "prediction" indices
% are along the columns. Element (i, j) shows in what fraction of cases in
% which the truth is j, the prediction was i.

M = [mean(preds(truths == 0) == 0) mean(preds(truths == 0) == 1) ;
     mean(preds(truths == 1) == 0) mean(preds(truths == 1) == 1)];

end