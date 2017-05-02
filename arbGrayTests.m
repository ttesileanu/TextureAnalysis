%% Check that binary stats are 10d

% make this reproducible
rng('default');

% patch settings
nPatches = 1000;
patchSize = 128;

% binary patches with biased entries (70% are 1)
binPatches = 1.0*(rand([patchSize, patchSize, nPatches]) < 0.7);

% stats for all patches
tic;
stats = [];
for i = 1:nPatches
    crtP = processBlock(binPatches(:, :, i));
    if isempty(stats)
        stats = zeros(length(crtP), nPatches);
    end
    stats(:, i) = crtP; %#ok<SAGROW>
end
disp(['Binary calculation took ' num2str(1000*toc, '%.0f'), ' ms.']);

% check dimensionality
stats0 = bsxfun(@minus, stats, mean(stats, 2));
[statsU, statsS, statsV] = svd(stats0, 0);

% there are 6 constraints (so the dimensionality should be 10), but they
% are not all exactly satisfied -- this is because the patch size is
% finite, so there are edge effects
figure;
plot(diag(statsS), '.k');
hold on;
plot([0 16], repmat(2.0/patchSize, 1, 2), '--k');
text(13, 2.0/patchSize, '2/patchSize', 'color', 'r', 'verticalalignment', 'bottom');
title('Principal component analysis (binary)');
xlabel('PC index');
ylabel('Singular value');

beautifygraph;
preparegraph;

% check the specific constraints that we expect to hold
[~, cons] = getStatistics;
values = zeros(size(cons, 1), nPatches);
for i = 1:size(cons, 1)
    values(i, :) = sum(bsxfun(@times, cons(i, :)', stats), 1);
end
expected = [0 0 0 0 0 1]';

figure;
plot([0 7], [0 0], ':k');
hold on;
plot(1:6, mean(values, 2) - expected, '.k');
errorbar(1:6, mean(values, 2) - expected, std(values, [], 2), 'marker', 'none', 'linestyle', 'none');

plot([0 7], repmat(1.0/patchSize, 1, 2), '--k');
plot([0 7], -repmat(1.0/patchSize, 1, 2), '--k');
text(5, 1.0/patchSize, '1/patchSize', 'color', 'r', 'verticalalignment', 'bottom');
text(5, -1.0/patchSize, '-1/patchSize', 'color', 'r', 'verticalalignment', 'top');
xlim([0 7]);
ylim([-1.5/patchSize, 1.5/patchSize]);
title('Constraints (binary)');
xlabel('Constraint index');
ylabel('Constraint value');

beautifygraph;
preparegraph;

%% Check that inf and binary stats agree on binary images

tic;
statsAgain = zeros(size(stats));
for i = 1:nPatches
    statsAgain(:, i) = processBlock(binPatches(:, :, i), inf);
end
disp(['Grayscale calculation took ' num2str(1000*toc, '%.0f'), ' ms.']);

disp(['shouldBeZero = ' num2str(norm(flatten(statsAgain - stats))) '.']);

%% What's the dimensionality on grayscale images?

grayPatches = abs(randn([patchSize, patchSize, nPatches]));

% stats for all patches
tic;
grayStats = [];
for i = 1:nPatches
    crtP = processBlock(grayPatches(:, :, i), inf);
    if isempty(grayStats)
        grayStats = zeros(length(crtP), nPatches);
    end
    grayStats(:, i) = crtP; %#ok<SAGROW>
end
disp(['Calculation on grayscale patches took ' num2str(1000*toc, '%.0f'), ' ms.']);

% check dimensionality
grayStats0 = bsxfun(@minus, grayStats, mean(grayStats, 2));
[grayStatsU, grayStatsS, grayStatsV] = svd(grayStats0, 0);

% there are 6 constraints (so the dimensionality should be 10), but they
% are not all exactly satisfied -- this is because the patch size is
% finite, so there are edge effects
figure;
plot(diag(grayStatsS), '.k');
hold on;
plot([0 16], repmat(2.0/patchSize, 1, 2), '--k');
text(13, 2.0/patchSize, '2/patchSize', 'color', 'r', 'verticalalignment', 'bottom');
title('Principal component analysis (grayscale)');
xlabel('PC index');
ylabel('Singular value');

beautifygraph;
preparegraph;

% check the specific constraints that used to hold in the binary case
[~, cons] = getStatistics;
grayValues = zeros(size(cons, 1), nPatches);
for i = 1:size(cons, 1)
    grayValues(i, :) = sum(bsxfun(@times, cons(i, :)', grayStats), 1);
end
expected = [0 0 0 0 0 1]';

figure;
plot([0 7], [0 0], ':k');
hold on;
plot(1:6, mean(grayValues, 2) - expected, '.k');
errorbar(1:6, mean(grayValues, 2) - expected, std(grayValues, [], 2), 'marker', 'none', 'linestyle', 'none');

plot([0 7], repmat(1.0/patchSize, 1, 2), '--k');
plot([0 7], -repmat(1.0/patchSize, 1, 2), '--k');
text(5, 1.0/patchSize, '1/patchSize', 'color', 'r', 'verticalalignment', 'bottom');
text(5, -1.0/patchSize, '-1/patchSize', 'color', 'r', 'verticalalignment', 'top');
xlim([0 7]);
ylim([-1.5/patchSize, 1.5/patchSize]);
title('Constraints (grayscale)');
xlabel('Constraint index');
ylabel('Constraint value');

beautifygraph;
preparegraph;