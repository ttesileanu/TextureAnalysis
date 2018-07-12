% generate the threshold predictions from the ternary natural image texture
% analysis

%% Load the ternary NI distribution

niAll = open(fullfile('save', 'TernaryDistribution_PennNoSky.mat'));

% choose one of the analyses
NRselection = [2, 32];
idx = find(cellfun(@(nr) isequal(nr, NRselection), niAll.valuesNR));
if isempty(idx)
    error('Can''t find NR selection in NI distribution file.');
end
if length(idx) > 1
    error('FOund multiple matches to the NR selection.');
end

ni0 = niAll.results{idx};

%% NI distribution preprocessing

% we typically restrict our analysis to patches identified as in-focus
restrictToFocus = true;

% check that the distribution we loaded has focus information
ni = rmfield(ni0, 'focus');
if restrictToFocus && isfield(ni, 'focus')
    disp('Restricting to in-focus patches.');
    mask = (ni0.focus.clusterIds == ni0.focus.focusCluster);
    fields = {'ev', 'patchLocations', 'imageIds'};
    for i = 1:length(fields)
        ni.(fields{i}) = ni.(fields{i})(mask, :);
    end
    ni.covM = cov(ni.ev);
end

%% Load the psychophysics data

ternaryAvg = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));
ternaryBySubject = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
    'subjects', '*', 'keepNaN', false);

%% Calculate gains and predicted thresholds

% use only second-order groups to set the overall scaling of the predictions
[gain, predictions, predictionDetails] = getPredictionsFromTernaryStats(...
    ni.ev, ternaryAvg, 'fitScaleOptions', {'mask', cellfun(@length, ternaryAvg.groups) == 6});

%% Save

save(fullfile('save', 'TernaryNIPredictions_PennNoSky_2x32.mat'), 'NRselection', ...
    'restrictToFocus', 'gain', 'predictions', 'predictionDetails');

%% Check match in single planes

plotTernaryMatrix(predictions, ternaryAvg, 'ellipses', false, 'mask', ~strcmp(ternaryAvg.groups, 'A_1'));

%% Check match in mixed planes

plotTernaryMatrix(predictions, ternaryAvg, 'ellipses', true, 'multi', 2);