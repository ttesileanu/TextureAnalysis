%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   dbChoice
%       Choice of database to use. Available options are
%           'PennNoSky' -- Penn image database, filtering out pictures with
%                          lots of sky
%           'vanHateren'    -- van Hateren image database, filtering out
%                              pictures with lots of sky or lots of
%                              human-made objects
%   compressType
%       Choose the way in which image values were compressed in the [0, 1]
%       interval before ternarizing. Options are
%           'equalize' -- histogram equalization
%           'contrast' -- contrast adaptation
%   NRselection
%       Choose one of the analyses, based on block-averaging factor (N) and
%       patch size (R).
%   symmetrizePP
%       Set to `true` to take the average between each psychophysics
%       measurement and the measurement in the opposite texture direction.
%       This effectively forces the measurements to be centered at the
%       origin.
%   gainTransform
%       A function to apply to the gains obtained from efficient coding.
%       This can be either a function handle or one of
%        'identity'
%           The gains are kept as they are.
%        'square'
%           The gains are squared. This was used in Hermundstad et al.,
%           leading to threshold predictions that are inversely proportional
%           to natural image standard deviations instead of their square
%           roots. Since the efficient coding problem solved here uses a
%           Gaussian approximation, this transformation might indicate a
%           departure of visual processing in the brain from Gaussianity.
%   nSamples
%       Number of samples to use for statistical tests.

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('NRselection', [2, 32]);
setdefault('symmetrizePP', false);
setdefault('gainTransform', 'square');
setdefault('nSamples', 10000);

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
% niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];
NRstr = [int2str(NRselection(1)) 'x' int2str(NRselection(2))];
niPredFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr ...
    '_' gainTransform '.mat'];

%% Load data

pp0 = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

if symmetrizePP
    % make sure data is symmetric
    ppOriginal = pp0;
    
    reflectTrafo = @(group) applyGroupReflection(group, 3);
    ppReflected = applyToThresholds(pp0, reflectTrafo, 'closed', true);
    
    pp0 = averageMeasurements(ppOriginal, ppReflected);
end

%% Load natural image predictions

niStructure = open(fullfile('save', niPredFileName));
ni0 = niStructure.predictions;

%% Select only second-order planes with finite thresholds

maskPP = cellfun(@(s) length(s) == 6 || sum(s == ';') == 1, pp0.groups) & ...
    isfinite(pp0.thresholds);

pp = selectMeasurements(pp0, maskPP);
ni = selectMeasurements(ni0, cellfun(@(s) length(s) == 6 || sum(s == ';') == 1, ni0.groups));

%% NULL MODEL #1: shuffle all thresholds randomly

% Note that this will in general not lead to elliptical threshold contours
% in texture planes.

% generate the samples
ppArtificialSimple = cell(1, nSamples);

comparisonsToExperimentSimple = zeros(nSamples, 1);
comparisonsToTheorySimple = zeros(nSamples, 1);

compType = 'direct';

progress = TextProgress;
for i = 1:nSamples
    ppArtificialSimple{i} = pp;
    permutation = randperm(length(ppArtificialSimple{i}.thresholds));
    ppArtificialSimple{i}.thresholds = ppArtificialSimple{i}.thresholds(permutation);
    
    % get distances to experiment and NI predictions, normalizing out the
    % median log threshold
    comparisonsToExperimentSimple(i) = compareMeasurements(pp, ppArtificialSimple{i}, compType, ...
        'normalize', true);
    comparisonsToTheorySimple(i) = compareMeasurements(ni, ppArtificialSimple{i}, compType, ...
        'normalize', true);
    
    progress.update(i*100/nSamples);
end
progress.done;

% distance between prediction and experiment
actualPPNIComparison = compareMeasurements(pp, ni, compType, 'normalize', true);

%% Make plots

fig = figure;
fig.Units = 'inches';
fig.Position = [4 3 7 5];

lo = min([comparisonsToExperimentSimple ; comparisonsToTheorySimple]);
hi = max([comparisonsToExperimentSimple ; comparisonsToTheorySimple]);
bins = linspace(lo, hi, 25);

subplot(2, 1, 1);
histogram(comparisonsToExperimentSimple, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to experiment');
beautifygraph;

subplot(2, 1, 2);
histogram(comparisonsToTheorySimple, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to theory');
beautifygraph;

preparegraph;

pvalSimple = mean(comparisonsToTheorySimple <= actualPPNIComparison);
if pvalSimple == 0
    disp(['p < 1/' int2str(nSamples) '.']);
else
    disp(['p = ' num2str(pvalSimple) '.']);
end
disp(['Actual D = ' num2str(actualPPNIComparison, '%.3f')]);
[lo, hi] = getHdi(comparisonsToTheorySimple, 0.95);
disp(['95% CI of D is [' num2str(lo, '%.3f') ', ' num2str(hi, '%.3f') ']']);

%% NULL MODEL #2: thresholds along 99 texture directions drawn randomly i.i.d.

% ...from log-normal distribution with mean and standard deviation taken
% from actual data. No correlations between directions.

% This naturally leads to elliptical threshold contours in every texture
% plane. At least one of the ellipse axes is always aligned with a
% coordinate axis.

% generate the samples
ppArtificial = cell(1, nSamples);
axisParams = [nanmean(log(pp.thresholds)), nanstd(log(pp.thresholds))];

comparisonsToExperiment = zeros(nSamples, 1);
comparisonsToTheory = zeros(nSamples, 1);

compType = 'direct';

progress = TextProgress;
for i = 1:nSamples
    ppArtificial{i} = inventMeasurements(pp.groups, pp.directions, 'axisParams', axisParams);
    
    % get distances to experiment and NI predictions, normalizing out the
    % median log threshold
    comparisonsToExperiment(i) = compareMeasurements(pp, ppArtificial{i}, compType, ...
        'normalize', true);
    comparisonsToTheory(i) = compareMeasurements(ni, ppArtificial{i}, compType, ...
        'normalize', true);
    
    progress.update(i*100/nSamples);
end
progress.done;

% distance between prediction and experiment
actualPPNIComparison = compareMeasurements(pp, ni, compType, 'normalize', true);

%% Make plots

fig = figure;
fig.Units = 'inches';
fig.Position = [4 3 7 5];

lo = min([comparisonsToExperiment ; comparisonsToTheory]);
hi = max([comparisonsToExperiment ; comparisonsToTheory]);
bins = linspace(lo, hi, 25);

subplot(2, 1, 1);
histogram(comparisonsToExperiment, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to experiment');
beautifygraph;

subplot(2, 1, 2);
histogram(comparisonsToTheory, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to theory');
beautifygraph;

preparegraph;

pval = mean(comparisonsToTheory <= actualPPNIComparison);
if pval == 0
    disp(['p < 1/' int2str(nSamples) '.']);
else
    disp(['p = ' num2str(pval) '.']);
end
disp(['Actual D = ' num2str(actualPPNIComparison, '%.3f')]);
[lo, hi] = getHdi(comparisonsToTheory, 0.95);
disp(['95% CI of D is [' num2str(lo, '%.3f') ', ' num2str(hi, '%.3f') ']']);

%% BAYESIAN MODEL: interpolate between constant threshold and NI predictions

% Consider a probabilistic model in which the log threshold in direction i
% is related to the natural image prediction in that direction by
%   log_th(i) - avg_log_th = beta*(log_pred - avg_log_pred) + eta(i)
% where eta(i) are i.i.d. errors drawn from a normal distribution with
% standard deviation sigma.

% Note that the log thresholds are simply random when beta = 0, and are
% given by the NI predictions + random error when beta = 1. We aim to find
% the posterior distribution for beta and see whether it overlaps 0 or 1.

% We try a few different priors to make sure we're not biasing the
% inference.

% prepare the data and NI predictions
% comparing mean-centered log thresholds
idxShuffle = matchMeasurements(pp, ni);
xi = log(pp.thresholds);
yi = log(ni.thresholds(idxShuffle));

xiHat = xi - mean(xi);
yiHat = yi - mean(yi);

% prepare the priors
n = length(xiHat);
priorBetaStd = 10;
priorLogSigmaStd = 10;
logPriors = {...
    {'flat', @(beta, logSigma) 0}, ...
    {'normalBeta', @(beta, logSigma) -0.5*(beta/priorBetaStd)^2}, ...
    {'normalLogSigma', @(beta, logSigma) -0.5*(logSigma/priorLogSigmaStd)^2}, ...
    {'normalBoth', @(beta, logSigma) -0.5*(beta/priorBetaStd)^2 - 0.5*(logSigma/priorLogSigmaStd)^2}};

% run the inference, plot some diagnostics
fig = figure;
posteriorSamplesBurned = cell(length(logPriors), 1);
% nSamples = 10000;
for i = 1:length(logPriors)
    priorName = logPriors{i}{1};
    logPrior = logPriors{i}{2};
    logPdf = @(params) -0.5*sum((params(1)*yiHat - xiHat).^2)/exp(2*params(2)) - ...
        n/2*log(2*pi) - n*params(2) + logPrior(params(1), params(2));
    
    % choose a neutral initial point
    initial = [0.5 0.0];
    posteriorSamples = slicesample(initial, nSamples, 'logpdf', logPdf);
    
    % plot diagnostics
    subplot(length(logPriors), 1, i);
    plot(posteriorSamples);
    legend({'\beta', 'log(\sigma)'});
    xlabel('Iteration');
    title(priorName);
    
    % store the samples
    posteriorSamplesBurned{i} = posteriorSamples(end/2:end, :);
end

% show posterior beta
minBeta = min([0 ; cellfun(@(samples) min(samples(:, 1)), posteriorSamplesBurned(:))]);
maxBeta = max([1 ; cellfun(@(samples) max(samples(:, 1)), posteriorSamplesBurned(:))]);
plotter = MatrixPlotter(length(logPriors), 'fixedSize', [8 4], 'fixedSizeUnits', 'inches');
while plotter.next
    i = plotter.index;
    
    histogram(posteriorSamplesBurned{i}(:, 1), 100, 'BinLimits', [minBeta maxBeta]);
    xlabel('\beta');
    title(logPriors{i}{1});
    
    [lo, hi] = getHdi(posteriorSamplesBurned{i}(:, 1), 0.95);
    disp(['95% credible interval for beta, ', logPriors{i}{1} ' prior: [' ...
        num2str(lo, '%.3f') ', ' num2str(hi, '%.3f') ']']);
end

% show posterior sigma
plotter = MatrixPlotter(length(logPriors), 'fixedSize', [8 4], 'fixedSizeUnits', 'inches');
while plotter.next
    i = plotter.index;
    
    histogram(exp(posteriorSamplesBurned{i}(:, 2)), 100, 'BinLimits', [0 0.5]);
    xlabel('\sigma');
    title(logPriors{i}{1});
    
    [lo, hi] = getHdi(exp(posteriorSamplesBurned{i}(:, 2)), 0.95);
    disp(['95% credible interval for sigma, ', logPriors{i}{1} ' prior: [' ...
        num2str(lo, '%.3f') ', ' num2str(hi, '%.3f') ']']);
end
