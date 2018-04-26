function plotBinaryThresholds(pred, meas, varargin)
% plotBinaryThresholds Make plot comparing predicted to measured thresholds
% for binary texture psychophysics.
%   plotBinaryThresholds(pred, meas) makes a plot comparing the predicted,
%   `pred`, to measured, `meas`, thresholds for binary texture
%   psychophysics. `meas` should be a structure as returned by
%   loadBinaryPP, containing cell arrays `groups` and `thresholds`. If the
%   structure also contains a `subjects` field, the plot will display
%   different subjects in different colors. `pred` can be a vector or a
%   cell array of vectors giving the predicted thresholds with the same
%   ordering as that used in `meas`. When there are several prediction
%   vectors available, they will be plotted in different colors. The
%   predictions can also be absent (set `pred` to an empty matrix).
%
%   The function will ignore multi-group entries.
%
%   Options:
%    'exclude'
%       Cell array of directions to exclude.
%    'collapsedirs'
%       Set to true to show only one measurement for both cardinal betas,
%       one for both diagonal betas, and one for all theta coordinates.
%    'symbolsize'
%       Size to use for the plot symbols.
%    'predlabels'
%       Labels to show for the predictions.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('exclude', {'A_1'}, @(c) iscell(c) && (isvector(c) || isempty(c)));
parser.addParameter('collapsedirs', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('symbolsize', 7, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('predlabels', {}, @(c) iscell(c) && (isvector(c) || isempty(c)));

if nargin == 1 && strcmp(pred, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% turn predictions into a cell array
if ~isempty(pred) && ~iscell(pred)
    pred = {pred};
end

% these are the groups we're looking at
edirs = {'A_1', 'AC_1_1', 'AB_1_1', 'AD_1_1', 'BC_1_1', ...
    'ABC_1_1_1', 'BCD_1_1_1', 'ABD_1_1_1', 'ACD_1_1_1', 'ABCD_1_1_1_1'};
% edirs = setdiff(edirs, params.exclude);

% read off sensitivities from predictions and measurements

% set up subjects
if isfield(meas, 'subjects')
    subjects = unique(meas.subjects);
    n_subjects = length(subjects);
else
    subjects = {'subj'};
    n_subjects = 1;
end
pred_mat = nan(length(pred), length(edirs));
meas_mat = nan(n_subjects, length(edirs));

for i = 1:length(edirs)
    % skip excluded directions
    if ismember(edirs{i}, params.exclude)
        continue;
    end
    
    % find which measurements correspond to this direction
    i0 = find(strcmp(meas.groups, edirs{i}));
    % nothing to do if none are found
    if isempty(i0)
        continue;
    end
    
    % iterate through the found thresholds
    for j = 1:length(i0)
        % assign thresholds to the appropriate subjects
        if isfield(meas, 'subjects')
            crt_subject = meas.subjects(i0(j));
        else
            crt_subject = subjects{1};
        end
        j0 = find(strcmp(subjects, crt_subject), 1);
        if ~isempty(j0)
            meas_mat(j0, i) = meas.thresholds(i0(j));
        end
    end
    
    % copy over the appropriate predictions
    for j0 = 1:length(pred)
        % keep only the first prediction for each group
        pred_mat(j0, i) = pred{j0}(i0(1));
    end
end

% sensitivities are inverse thresholds
pred_sens = 1./pred_mat;
meas_sens = 1./meas_mat;

% labels for the coordinates
c_labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% figure out where to plot things -- where the labels are, where the
% predictions are, where the measurements are
% XXX this is cumbersome...
if params.collapsedirs
    if ~isempty(pred)
        if ismember('A_1', params.exclude)
            label_positions = [1:2 5:6 9:12 15];
            meas_positions = [3 7 13 16];
        else
            label_positions = [-2 1:2 5:6 9:12 15];
            meas_positions = [-1 3 7 13 16];
        end
        pred_positions = label_positions;
        label_choices = c_labels((1 + ismember('A_1', params.exclude)):end);
    else
        if ismember('A_1', params.exclude)
            meas_positions = [1 4 7 10];
            label_choices = {};
        else
            meas_positions = [-2 1 4 7 10];
            label_choices = c_labels(1);
        end
        pred_positions = [];
        label_positions = meas_positions;
        label_choices = [label_choices ...
            {[c_labels{2} ',' c_labels{3}] ...
            [c_labels{4} ',' c_labels{5}] ...
            [c_labels{6} ',' c_labels{7} ',' c_labels{8} ',' c_labels{9}] ...
            c_labels{10}}];
    end
else
    if ismember('A_1', params.exclude)
        label_positions = 1:3:27;
    else
        label_positions = -2:3:27;
    end
    if ~isempty(pred)
        pred_positions = label_positions - 1;
        meas_positions = label_positions + 1;
    else
        pred_positions = [];
        meas_positions = label_positions;
    end
    label_choices = c_labels((1 + ismember('A_1', params.exclude)):end);
end

was_hold = ishold;
hold on;

% draw the predictions
if ~isempty(pred)
    pred_colors = cool(size(pred_sens, 1));
    for i = 1:size(pred_sens, 1)
        % skip luminance if necessary
        if ismember('A_1', params.exclude)
            crt_pred_sens = pred_sens(i, 2:end);
        else
            crt_pred_sens = pred_sens(i, :);
        end
        
        % draw
        plot(pred_positions, crt_pred_sens, 'o', ...
            'markersize', params.symbolsize, 'color', pred_colors(i, :), 'linewidth', 2);
    end
end

% draw the measurements
if size(meas_sens, 1) > 1
    select_half = @(m) m(1:end/2, :);
    pp_colors = select_half(hot(2*size(meas_sens, 1)));
else
    pp_colors = [1 0 0];
end
for i = 1:size(meas_sens, 1)
    crt_meas_sens0 = meas_sens(i, :);
    
    % collapse if necessary
    if params.collapsedirs
        crt_meas_sens = [crt_meas_sens0(1) geomean(crt_meas_sens0(2:3)) ...
            geomean(crt_meas_sens0(4:5)) geomean(crt_meas_sens0(6:9)) ...
            crt_meas_sens0(10)];
    else
        crt_meas_sens = crt_meas_sens0;
    end
    % skip luminance if necessary
    if ismember('A_1', params.exclude)
        crt_meas_sens = crt_meas_sens(2:end);
    end
    
    % draw 
    plot(meas_positions, crt_meas_sens, ...
        's', 'color', pp_colors(i, :), 'markersize', params.symbolsize, 'linewidth', 2);
end


% set up labels and ticks
title('Texture sensitivity from natural images vs. psychophysics');
ylabel('Sensitivity');

xlim([min(label_positions) - 1, max(label_positions) + 2]);
ylim([0 ceil(2*nanmax([pred_sens(:) ; meas_sens(:)]))/2]);

% set up legend
if length(subjects) > 1
    pp_labels = cellfun(@(s) ['data ' s], subjects, 'uniform', false);
else
    pp_labels = {'data'};
end

if isempty(params.predlabels)
    ni_labels = arrayfun(@(i) ['pred ' int2str(i)], 1:size(pred_sens, 1), 'uniform', false);
else
    ni_labels = params.predlabels;
end

legend([ni_labels(:) ; pp_labels(:)]);

beautifygraph;

set(gca, 'xtick', label_positions, 'xticklabel', label_choices, 'xminortick', 'off');

if ~was_hold
    hold off;
end

end