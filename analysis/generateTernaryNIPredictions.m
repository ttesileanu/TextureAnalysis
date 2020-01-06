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
%   filterScope
%       Choose the scope of the whitening filter:
%           'patch' (default) -- whiten each patch separately
%           'image'           -- whiten before patchifying
%   compressScope
%       Choose the scope of the compression. Note that setting this to
%       'image' forces `filterScope` to be 'image' as well.
%           'patch' (default) -- compress each patch separately
%           'image'           -- compress before patchifying
%   NRselection
%       Choose one of the analyses, based on block-averaging factor (N) and
%       patch size (R). This can also be a cell array of tuples in order to
%       generate several analyses at the same time.
%   restrictToFocus
%       Set to `true` to only keep patches that were identified as in-focus
%       by a two-Gaussian fit.
%   useErrorBars
%       Set to `false` to ignore error bars when finding the optimal
%       scaling for NI predictions.
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
%       If `true`, use a second fitting parameter to estimate the slope
%       between log predictions and log measurements.

setdefault('dbChoice', 'PennNoSky');
setdefault('compressType', 'equalize');
setdefault('compressScope', 'patch');
setdefault('filterScope', 'patch');
setdefault('NRselection', {[1, 32], [1, 48], [1, 64], [2, 32], [2, 48], [2, 64], ...
    [4, 32], [4, 48], [4, 64]});
setdefault('restrictToFocus', true);
setdefault('useErrorBars', true);
setdefault('gainTransform', 'square');
setdefault('fitLogSlope', false);

%% Preprocess options

if strcmp(compressScope, 'image')
    filterScope = 'image';
end

if ~strcmp(compressType, 'equalize')
    extras = ['_' compressType];
else
    extras = '';
end
if ~strcmp(filterScope, 'patch')
    extras = [extras '_flt' filterScope];
end
if ~strcmp(compressScope, 'patch')
    extras = [extras '_comp' compressScope];
end
niFileName = ['TernaryDistribution_' dbChoice extras '.mat'];

switch gainTransform
    case 'identity'
        gainTransformFct = @(x) x;
    case 'square'
        gainTransformFct = @(x) x.^2;
    otherwise
        if ~isfa(gainTransform, 'function_handle')
            error('Invalid value for gainTransform.');
        end
        gainTransformFct = gainTransform;
end

% handle running multiple analyses at the same time
if iscell(NRselection)
    allNRs__ = NRselection;s
    for selIt__ = 1:length(allNRs__)
        NRselection = allNRs__{selIt__};
        disp(['Working on N=',s int2str(NRselection(1)), ', R=' ...
            int2str(NRselection(2)) '...']);
        clearvars -except dbChoice compressType compressScope filterScope ...
            NRselection restrictToFocus useErrorBars gainTransform fitLogSlope ...
            allNRs__ selIt__;
        generateTernaryNIPredictions;
    end
    return;
end

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
% pp_extra = open('data/extra_ternary_thresholds.mat');
% pp_extra_AC12 = selectMeasurements(pp_extra.avg, ...
%     strcmp(pp_extra.avg.groups, 'AC_1_2'));
% 
% pp = catMeasurements(pp, pp_extra_AC12);

%% Calculate gains and predicted thresholds

% use only second-order groups to set the overall scaling of the predictions
if useErrorBars
    ppForFit = pp;
else
    ppForFit = rmfield(pp, 'thresholdIntervals');
end
% [gain, predictions, predictionDetails] = getPredictionsFromTernaryStats(...
%     ni.ev, ppForFit, 'fitScaleOptions', {'mask', cellfun(@length, pp.groups) == 6});
[gain, predictions, predictionDetails] = getPredictionsFromTernaryStats(...
    ni.ev, ppForFit, ...
    'fitScaleOptions', {'mask', cellfun(@(s) length(s) == 6 || sum(s == ';') == 1, pp.groups), ...
        'logSlope', fitLogSlope}, ...
    'efficientCodingOptions', {'gainTransform', gainTransformFct});

%% Save

NRstr = [int2str(NRselection(1)) 'x' int2str(NRselection(2))];
fitLogSuffixes = {'', '_powfit'};
outFileName = ['TernaryNIPredictions_' dbChoice extras '_' NRstr ...
    '_' gainTransform fitLogSuffixes{1 + fitLogSlope} '.mat'];
save(fullfile('save', outFileName), 'NRselection', 'restrictToFocus', ...
    'gain', 'predictions', 'predictionDetails', 'gainTransformFct');

%% Check match

fig = figure;
fig.Units = 'inches';
totalX = 12;
totalY = 10;
fig.Position = [2 2 totalX totalY];

% axes for mixed planes
multi_ax = zeros(22, 1);
rowNumbers = [2 2 6 6 6];
rowShifts = [4 4 0 0 0];
rowStarts = [1 1+cumsum(rowNumbers(1:end-1))];
figX = totalX/max(rowNumbers);
figY = totalY/length(rowNumbers);
factorX = 0.75;
factorY = 0.75;
for i = 1:length(multi_ax)
    crtAx = axes;
    
    crtRow = find(i >= rowStarts, 1, 'last');
    crtShift = rowShifts(crtRow);
    crtCol = i - rowStarts(crtRow); % 0-based!!
    
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [(crtShift + crtCol)*figX totalY - crtRow*figY figX*factorX figY*factorY];
    
    multi_ax(i) = crtAx;
end

% axes for single planes
single_ax = zeros(4, 1);
for i = 1:length(single_ax)
    crtAx = axes;
    
    crtRow = floor((i-1)/2);
    crtCol = mod(i-1, 2);
    
    crtAx.Units = 'inches';
    crtAx.OuterPosition = [crtCol*2*figX totalY - (1 + crtRow)*figY figX figY];
    
    single_ax(i) = crtAx;
end

% plot single planes
plusMinus = '+-';
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) length(g) == 6 && ~strcmp(g(1:2), 'AC'), ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'triangleOptions', {'fontscale', 0.667, 'edgelabels', 'digit'}, ...
    'limits', 1.5, ...
    'plotterOptions', {'fixedAxes', single_ax}, ...
    'titleShift', [0 -0.5], 'titleAlignment', {'center', 'bottom'}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)});

% plot mixed planes
plotTernaryMatrix({predictions, pp}, 'ellipse', false, ...
    'groupMaskFct', @(g) sum(g == ';') == 1, ...
    'beautifyOptions', {'ticks', 'off', 'ticklabels', false, ...
        'titlesize', 12, 'titleweight', 'normal', 'noaxes', true, ...
        'fontscale', 0.667}, ...
    'plotterOptions', {'fixedAxes', multi_ax}, ...
    'labelOptions', {'FontSize', 8, 'subscriptSpacing', -0.55, ...
        'coeffToStr', @(i) plusMinus(i)}, ...
    'xLabelAlignment', {'center', 'bottom'}, 'yLabelAlignment', {'left', 'bottom'});
preparegraph;
