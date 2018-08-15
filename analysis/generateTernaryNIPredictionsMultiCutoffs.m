% calculate predictions and prediction accuracy as a function of
% ternarization cut points

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   restrictToFocus
%       Set to `true` to only keep patches that were identified as in-focus
%       by a two-Gaussian fit.

setdefault('restrictToFocus', true);

%% Load the psychophysics data

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));

% add additional data from Jonathan, but keep only AC_1_2 plane
pp_extra = open('data/extra_ternary_thresholds.mat');
pp_extra_AC12 = selectMeasurements(pp_extra.avg, ...
    strcmp(pp_extra.avg.groups, 'AC_1_2'));

pp = catMeasurements(pp, pp_extra_AC12);

%% Load ternary distributions, preprocess, and get predictions

allPredictionResults = struct;

allFiles = dir(fullfile('save', 'multicutoff'));

crt = 1;
for i = 1:length(allFiles)
    crtFile = allFiles(i);
    
    % keep only files (not folders) ending in .mat
    if crtFile.isdir
        continue;
    end
    
    % XXX should restrict to a common prefix
    crtFile = fullfile(crtFile.folder, crtFile.name);
    [~, ~, crtExt] = fileparts(crtFile);
    if ~strcmp(crtExt, '.mat')
        continue;
    end
    
    disp(['Processing ' crtFile '...']);
    
    % load data, restrict to in-focus patches if necessary
    crtData = open(crtFile);
    crtCutoffs = crtData.cutoffs;
    if isfield(crtData.results, 'focus')
        crtResults = rmfield(crtData.results, 'focus');
    else
        crtResults = crtData.results;
    end
    if restrictToFocus
        if ~isfield(crtData.results, 'focus') || ~isfield(crtData.results.focus, 'clusterIds')
            warning('%s does not contain focus information. Skipping.', crtFile);
            continue;
        end
        mask = (crtData.results.focus.clusterIds == crtData.results.focus.focusCluster);
        fields = {'ev', 'patchLocations', 'imageIds'};
        for j = 1:length(fields)
            crtResults.(fields{j}) = crtResults.(fields{j})(mask, :);
        end
        crtResults.covM = cov(crtResults.ev);
    end
    
    % use only second-order groups to set the overall scaling of the predictions
    [crtGain, crtPredictions, crtPredictionDetails] = getPredictionsFromTernaryStats(...
        crtResults.ev, pp, 'fitScaleOptions', {'mask', cellfun(@length, pp.groups) == 6});
    
    % store the results
    allPredictionResults(crt).file = crtFile;
    allPredictionResults(crt).cutoffs = crtCutoffs;
    allPredictionResults(crt).gain = crtGain;
    allPredictionResults(crt).predictions = crtPredictions;
    allPredictionResults(crt).predictionDetails = crtPredictionDetails;
    crt = crt + 1;
end

% sort the prediction results in increaing order of cutoffs(1)
[~, reorder] = sort(arrayfun(@(s) s.cutoffs(1), allPredictionResults));
allPredictionResults = allPredictionResults(reorder);

%% Save the results

save(fullfile('save', 'TernaryDistribution_PennNoSky_2x32_multicutoff_predictions.mat'));

%% Check how match to psychophysics depends on cutoff

grayAmounts = arrayfun(@(s) diff(s.cutoffs), allPredictionResults);
% restrict comparison to either second-order or mixed groups
groupMaskFct = @(g) length(g) == 6 || sum(g == ';') == 1;
differences = arrayfun(@(s) compareMeasurements(s.predictions, pp, 'group', ...
    'groupMaskFct', groupMaskFct), allPredictionResults);

plot(grayAmounts, differences, '.-');
xlabel('Fraction of gray');
ylabel('Mismatch between NI and PP');