% check the robustness of the two-gaussian fit by resampling with
% replacement

%% Load the binary NI statistics

niBinary = open('save/BinaryDistribution_PennNoSky.mat');

%% Load resampled binary stats

load(fullfile('save', 'BinaryFocusResamples_PennNoSky.mat'));

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

%% Load resampled ternary stats

load(fullfile('save', 'TernaryFocusResamples_PennNoSky.mat'));

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
