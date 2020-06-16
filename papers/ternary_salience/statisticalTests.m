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
%       patch size (R). This can also be a cell array of tuples in order to
%       generate several analyses at the same time.
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
%   fitLogSlope
%       Whether to use the predictions where a second fitting parameter was
%       employed to estimate the slope between log predictions and log
%       measurements.
%   nSamples
%       Number of samples to use for statistical tests.
%   outFile
%       Text file to which to append the results. Format, e.g.,
%         '2x48, 0.231, 0.0000, [0.512, 0.670], 0.0003, [0.412, 0.732], [0.451, 0.632], [0.241, 0.283]
%       where the first entry identifies the preprocessing options, next we
%       have the D value for the actual predictions obtained from our
%       natural-image analysis; then the p-value for the null hypothesis
%       used in the first permutation test (see the paper) and the 95%
%       highest-density interval of D values obtained from the permutations;
%       then the same for the second permutation test (see the paper); and
%       finally the highest-density interval for the posterior probability
%       of the beta and sigma parameters from the parameter-estimation test
%       (see the paper). These are the values that appear in, e.g., Table 1
%       in the SI. Note that a p-value of 0 saved to file means that none
%       of the `nSamples` samples reached a value as extreme as the
%       observed one -- thus, it is more appropriately summarized as
%       `p < 1 / nSamples`.
%
%       By default, the file chosen is located in the 'save' folder, and is
%       named, e.g., 'PennNoSky_equalize_square.csv', identifying the
%       database, compression type for lightness values, and gain transform
%       for the efficient coding results. The file is created if it doesn't
%       exist, but otherwise it is only appended to.
%   errorMeasure
%       Kind of measure to use for estimating errors. This can be
%        'RMSLE':   root mean square log error
%        'MedALE':  median absolute log error

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('NRselection', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});
setdefault('symmetrizePP', false);
setdefault('gainTransform', 'square');
setdefault('nSamples', 10000);
setdefault('fitLogSlope', false);
if ischar(gainTransform)
    gainTransformStr = ['_' gainTransform];
else
    gainTransformStr = '';
end
setdefault('outFile', fullfile('save', [dbChoice '_' compressType gainTransformStr '.csv']));
setdefault('errorMeasure', 'RMSLE');

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end

% handle running multiple analyses at the same time
if iscell(NRselection)
    allNRs__ = NRselection;
    for selIt__ = 1:length(allNRs__)
        NRselection = allNRs__{selIt__};
        disp(['Working on N=', int2str(NRselection(1)), ', R=' ...
            int2str(NRselection(2)) '...']);
        clearvars -except dbChoice compressType NRselection ...
            symmetrizePP gainTransform nSamples fitLogSlope outFile ...
            errorMeasure allNRs__ selIt__;
        statisticalTests;
    end
    return;
end

% niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];
NRstr = [int2str(NRselection(1)) 'x' int2str(NRselection(2))];
fitLogSuffixes = {'', '_powfit'};
niPredFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr ...
    '_' gainTransform fitLogSuffixes{1 + fitLogSlope} '.mat'];

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

% ensure reproducibility
rng(2873834);

% generate the samples
ppArtificialSimple = cell(1, nSamples);

comparisonsToExperimentSimple = zeros(nSamples, 1);
comparisonsToTheorySimple = zeros(nSamples, 1);

compOpts = {'direct', 'normalize', true, 'summaryType', errorMeasure};

progress = TextProgress;
for i = 1:nSamples
    ppArtificialSimple{i} = shuffleMeasurements(pp);
    
    % get distances to experiment and NI predictions, normalizing out the
    % median log threshold
    comparisonsToExperimentSimple(i) = compareMeasurements(pp, ppArtificialSimple{i}, compOpts{:});
    comparisonsToTheorySimple(i) = compareMeasurements(ni, ppArtificialSimple{i}, compOpts{:});
    
    progress.update(i*100/nSamples);
end
progress.done;

% distance between prediction and experiment
actualPPNIComparison = compareMeasurements(pp, ni, compOpts{:});

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

D1 = actualPPNIComparison;
D1lo = lo;
D1hi = hi;
p1 = pvalSimple;

%% NULL MODEL #2: shuffle all groups randomly

% ensure reproducibility
rng(4736734);

% generate the samples
ppArtificialPermGroupCyclic = cell(1, nSamples);

comparisonsToExperimentPermGroupCyclic = zeros(nSamples, 1);
comparisonsToTheoryPermGroupCyclic = zeros(nSamples, 1);

compOpts = {'direct', 'normalize', true, 'summaryType', errorMeasure};

progress = TextProgress;
for i = 1:nSamples
    ppArtificialPermGroupCyclic{i} = shuffleMeasurements(pp, 'group', true, ...
        'cyclic', true);
    
    % get distances to experiment and NI predictions, normalizing out the
    % median log threshold
    comparisonsToExperimentPermGroupCyclic(i) = compareMeasurements(...
        pp, ppArtificialPermGroupCyclic{i}, compOpts{:});
    comparisonsToTheoryPermGroupCyclic(i) = compareMeasurements(...
        ni, ppArtificialPermGroupCyclic{i}, compOpts{:});
    
    progress.update(i*100/nSamples);
end
progress.done;

% distance between prediction and experiment
actualPPNIComparison = compareMeasurements(pp, ni, compOpts{:});

%% Make plots

fig = figure;
fig.Units = 'inches';
fig.Position = [4 3 7 5];

lo = min([comparisonsToExperimentPermGroupCyclic ; comparisonsToTheoryPermGroupCyclic]);
hi = max([comparisonsToExperimentPermGroupCyclic ; comparisonsToTheoryPermGroupCyclic]);
bins = linspace(lo, hi, 25);

subplot(2, 1, 1);
histogram(comparisonsToExperimentPermGroupCyclic, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to experiment');
beautifygraph;

subplot(2, 1, 2);
histogram(comparisonsToTheoryPermGroupCyclic, bins);
hold on;
plot([actualPPNIComparison actualPPNIComparison], ylim, 'r--');
xlabel('comparison to theory');
beautifygraph;

preparegraph;

pvalSimple = mean(comparisonsToTheoryPermGroupCyclic <= actualPPNIComparison);
if pvalSimple == 0
    disp(['p < 1/' int2str(nSamples) '.']);
else
    disp(['p = ' num2str(pvalSimple) '.']);
end
disp(['Actual D = ' num2str(actualPPNIComparison, '%.3f')]);
[lo, hi] = getHdi(comparisonsToTheoryPermGroupCyclic, 0.95);
disp(['95% CI of D is [' num2str(lo, '%.3f') ', ' num2str(hi, '%.3f') ']']);

D2 = actualPPNIComparison;
D2lo = lo;
D2hi = hi;
p2 = pvalSimple;

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

% ensure reproducibility
rng(75236);

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
figure;
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
    
    if strcmp(logPriors{i}{1}, 'flat')
        beta_lo = lo;
        beta_hi = hi;
    end
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
    if strcmp(logPriors{i}{1}, 'flat')
        sigma_lo = lo;
        sigma_hi = hi;
    end
end

%% Save to file, if a file is provided

if ~isempty(outFile)
    f = fopen(outFile, 'a');
    fprintf(f, '%dx%d, %.3f, %.4f, [%.3f, %.3f], %.4f, [%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f]\n', ...
        NRselection(1), NRselection(2), D1, ...
        p1, D1lo, D1hi, ....
        p2, D2lo, D2hi, ...
        beta_lo, beta_hi, sigma_lo, sigma_hi);
    fclose(f); 
end

%% Make a pretty plot of the beta posterior

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [4 2];

fig.Color = [1 1 1];

idx = 1;

% histogram(posteriorSamplesBurned{idx}(:, 1), 'normalization', 'pdf');

histogram(posteriorSamplesBurned{idx}(:, 1), 50, 'BinLimits', [0.6 1.2], ...
    'normalization', 'pdf');
xlabel('\beta');
ylabel('Posterior PDF');

beautifygraph('linewidth', 0.5, 'minorticks', 'off');
preparegraph;
