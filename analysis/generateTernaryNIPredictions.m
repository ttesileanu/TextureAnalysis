% generate the threshold predictions from the ternary natural image texture
% analysis

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   dbChoice
%       Choice of database to use. Available options are
%           'PennNoSky' -- Penn image database, filtering out pictures with
%                          lots of sky
%   compressType
%       Choose the way in which image values were compressed in the [0, 1]
%       interval before ternarizing. Options are
%           'equalize' -- histogram equalization
%           'contrast' -- contrast adaptation
%   NRselection
%       Choose one of the analyses, based on block-averaging factor (N) and
%       patch size (R).
%   restrictToFocus
%       Set to `true` to only keep patches that were identified as in-focus
%       by a two-Gaussian fit.

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('NRselection', [2, 32]);
setdefault('restrictToFocus', true);

%% Preprocess options

if ~strcmp(compressType, 'equalize')
    compressExt = ['_' compressType];
else
    compressExt = '';
end
niFileName = ['TernaryDistribution_' dbChoice compressExt '.mat'];

%% Load the ternary NI distribution

niAll = open(fullfile('save', niFileName));

% choose one of the analyses
idx = find(cellfun(@(nr) isequal(nr, NRselection), niAll.valuesNR));
if isempty(idx)
    error('Can''t find NR selection in NI distribution file.');
end
if length(idx) > 1
    error('Found multiple matches to the NR selection.');
end

ni0 = niAll.results{idx};

%% NI distribution preprocessing

% check that the distribution we loaded has focus information
ni = rmfield(ni0, 'focus');
if restrictToFocus && isfield(ni0, 'focus')
    disp('Restricting to in-focus patches.');
    mask = (ni0.focus.clusterIds == ni0.focus.focusCluster);
    fields = {'ev', 'patchLocations', 'imageIds'};
    for i = 1:length(fields)
        ni.(fields{i}) = ni.(fields{i})(mask, :);
    end
    ni.covM = cov(ni.ev);
end

%% Load the psychophysics data

pp = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'));
% ternaryBySubject = loadTernaryPP(fullfile('data', 'mtc_soid_xlsrun_summ.mat'), ...
%     'subjects', '*', 'keepNaN', false);

% add additional data from Jonathan, but keep only AC_1_2 plane
pp_extra = open('data/extra_ternary_thresholds.mat');
pp_extra_AC12 = selectMeasurements(pp_extra.avg, ...
    strcmp(pp_extra.avg.groups, 'AC_1_2'));

pp = catMeasurements(pp, pp_extra_AC12);

%% Calculate gains and predicted thresholds

% use only second-order groups to set the overall scaling of the predictions
[gain, predictions, predictionDetails] = getPredictionsFromTernaryStats(...
    ni.ev, pp, 'fitScaleOptions', {'mask', cellfun(@length, pp.groups) == 6});

%% Save

NRstr = [int2str(NRselection(1)) 'x' int2str(NRselection(2))];
outFileName = ['TernaryNIPredictions_' dbChoice compressExt '_' NRstr '.mat'];
save(fullfile('save', outFileName), 'NRselection', 'restrictToFocus', ...
    'gain', 'predictions', 'predictionDetails');

%% Check match in single planes

plotTernaryMatrix({predictions, pp}, ...
    'ellipse', false, ...
    'groupMaskFct', @(group) ~strcmp(group, 'A_1') && sum(group == ';') == 0);

%% Check match in mixed planes

plotTernaryMatrix({predictions, pp}, ...
    'ellipse', true, ...
    'groupMaskFct', @(group) sum(group == ';') == 1);
