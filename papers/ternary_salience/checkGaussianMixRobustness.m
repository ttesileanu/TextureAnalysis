% check the robustness of the two-gaussian fit by resampling with
% replacement

%% Load the binary NI statistics

niBinary = open('save/BinaryDistribution_PennNoSky.mat');

%% Resample binary stats

% keep this reproducible
rng(34788);

nBinarySamples = 6;
nBinaryResults = length(niBinary.results);

resampledBinaryClusterInfo = repmat(...
    {struct('clusterIds', cell(nBinarySamples, 1), 'clusterDistances', cell(nBinarySamples, 1), ...
     'shuffleIdxs', cell(nBinarySamples, 1))}, ...
    nBinaryResults, 1);

t0 = tic;
for k = 1:nBinarySamples
    disp(['Sample ' int2str(k) ' / ' int2str(nBinarySamples) '... (' ...
        num2str(toc(t0), '%.2f') ' seconds)']);
    for i = 1:nBinaryResults
        crtResults = niBinary.results{i};
        sizeResults = size(crtResults.ev, 1);
        % resample
        crtShuffleIdxs = randi(sizeResults, sizeResults, 1);
        crtResampled = sampleResults(crtResults, crtShuffleIdxs);
        % recalculate 2-gaussian fit, without redoing sharpness calculation
        crtResampledWithFocus = runFocusAnalysis(crtResampled, 'reuseSharpness', true);
        % exchange cluster Ids if necessary
        originalFocusCluster = crtResampled.focus.focusCluster;
        if crtResampledWithFocus.focus.focusCluster ~= originalFocusCluster
            crtResampledWithFocus.focus.clusterIds = ...
                3 - crtResampledWithFocus.focus.clusterIds;
            crtResampledWithFocus.focus.clusterDistances = ...
                fliplr(crtResampledWithFocus.focus.clusterDistances);
            % no need to exchange other members because we don't use them
        end
        % store the results
        resampledBinaryClusterInfo{i}(k).clusterIds = crtResampledWithFocus.focus.clusterIds;
        resampledBinaryClusterInfo{i}(k).clusterDistances = crtResampledWithFocus.focus.clusterDistances;
        resampledBinaryClusterInfo{i}(k).shuffleIdxs = crtShuffleIdxs;
    end
end
disp(['Took ' num2str(toc(t0), '%.2f') ' seconds.']);

%% Load resampled binary stats

load(fullfile('save', 'BinaryFocusResamples_PennNoSky.mat'));

%% Save resampled binary stats

save(fullfile('save', 'BinaryFocusResamples_PennNoSky.mat'), ...
    'nBinaryResults', 'nBinarySamples', 'resampledBinaryClusterInfo');

%% Check robustness binary

% similarities of cluster IDs
similarities = zeros(nBinaryResults, nBinarySamples);
% similarities of cluster distances
distanceSimilarities = zeros(nBinaryResults, nBinarySamples);
for i = 1:nBinaryResults
    % for each pair of down-sampling ratio and patch size, find how similar
    % the cluster assignments are between repetitions of the 2-gaussian fit
    for k = 1:nBinarySamples
        crtShuffleIdxs = resampledBinaryClusterInfo{i}(k).shuffleIdxs;
        resampledClusterIds = resampledBinaryClusterInfo{i}(k).clusterIds;
        shuffledBaseClusterIdxs = niBinary.results{i}.focus.clusterIds(crtShuffleIdxs);
        similarities(i, k) = mean(resampledClusterIds == shuffledBaseClusterIdxs);
        
        resampledDistances = resampledBinaryClusterInfo{i}(k).clusterDistances;
        shuffledBaseDistances = niBinary.results{i}.focus.clusterDistances(crtShuffleIdxs, :);
        distanceSimilarities(i, k) = corr(resampledDistances(:), shuffledBaseDistances(:));
    end
    crtN = niBinary.valuesN(i);
    crtR = niBinary.valuesR(i);
    fprintf('N=%d, R=%d clusterIds similarities: %s\n', crtN, crtR, ...
        num2str(similarities(i, :), '%8.4f'));
    fprintf('N=%d, R=%d clusterDistances similarities: %s\n', crtN, crtR, ...
        num2str(distanceSimilarities(i, :), '%8.4f'));
end

%% Load the ternary NI statistics

niTernary = open('save/TernaryDistribution_PennNoSky.mat');

%% Resample ternary stats

% keep this reproducible
rng(459485);

nSamples = 6;
nResults = length(niTernary.results);

resampledClusterInfo = repmat(...
    {struct('clusterIds', cell(nSamples, 1), 'clusterDistances', cell(nSamples, 1), ...
     'shuffleIdxs', cell(nSamples, 1))}, ...
    nResults, 1);

t0 = tic;
for k = 1:nSamples
    disp(['Sample ' int2str(k) ' / ' int2str(nSamples) '... (' ...
        num2str(toc(t0), '%.2f') ' seconds)']);
    for i = 1:nResults
        crtResults = niTernary.results{i};
        sizeResults = size(crtResults.ev, 1);
        % resample
        crtShuffleIdxs = randi(sizeResults, sizeResults, 1);
        crtResampled = sampleResults(crtResults, crtShuffleIdxs);
        % recalculate 2-gaussian fit, without redoing sharpness calculation
        crtResampledWithFocus = runFocusAnalysis(crtResampled, 'reuseSharpness', true);
        % exchange cluster Ids if necessary
        originalFocusCluster = crtResampled.focus.focusCluster;
        if crtResampledWithFocus.focus.focusCluster ~= originalFocusCluster
            crtResampledWithFocus.focus.clusterIds = ...
                3 - crtResampledWithFocus.focus.clusterIds;
            crtResampledWithFocus.focus.clusterDistances = ...
                fliplr(crtResampledWithFocus.focus.clusterDistances);
            % no need to exchange other members because we don't use them
        end
        % store the results
        resampledClusterInfo{i}(k).clusterIds = crtResampledWithFocus.focus.clusterIds;
        resampledClusterInfo{i}(k).clusterDistances = crtResampledWithFocus.focus.clusterDistances;
        resampledClusterInfo{i}(k).shuffleIdxs = crtShuffleIdxs;
    end
end
disp(['Took ' num2str(toc(t0), '%.2f') ' seconds.']);

%% Load resampled ternary stats

load(fullfile('save', 'TernaryFocusResamples_PennNoSky.mat'));

%% Save resampled ternary stats

save(fullfile('save', 'TernaryFocusResamples_PennNoSky.mat'), ...
    'nResults', 'nSamples', 'resampledClusterInfo');

%% Check robustness ternary

similarities = zeros(nResults, nSamples);
for i = 1:nResults
    % for each pair of down-sampling ratio and patch size, find how similar
    % the cluster assignments are between repetitions of the 2-gaussian fit
    for k = 1:nSamples
        crtShuffleIdxs = resampledClusterInfo{i}(k).shuffleIdxs;
        resampledClusterIds = resampledClusterInfo{i}(k).clusterIds;
        shuffledBaseClusterIdxs = niTernary.results{i}.focus.clusterIds(crtShuffleIdxs);
        similarities(i, k) = mean(resampledClusterIds == shuffledBaseClusterIdxs);
        
        resampledDistances = resampledClusterInfo{i}(k).clusterDistances;
        shuffledBaseDistances = niTernary.results{i}.focus.clusterDistances(crtShuffleIdxs, :);
        distanceSimilarities(i, k) = corr(resampledDistances(:), shuffledBaseDistances(:));
    end
    crtN = niTernary.valuesN(i);
    crtR = niTernary.valuesR(i);
    fprintf('N=%d, R=%d clusterIds similarities: %s\n', crtN, crtR, ...
        num2str(similarities(i, :), '%8.4f'));
    fprintf('N=%d, R=%d clusterDistances similarities: %s\n', crtN, crtR, ...
        num2str(distanceSimilarities(i, :), '%8.4f'));
end
