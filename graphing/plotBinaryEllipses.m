function plotBinaryEllipses(pred, meas, varargin)
% plotBinaryEllipses Make plot comparing predicted to measured threshold
% ellipses in mixed planes for binary texture psychophysics.
%   plotBinaryEllipses(pred, meas) makes a plot comparing the predicted,
%   `pred`, to measured, `meas`, threshold ellipses obtained in mixed
%   binary texture groups. `meas` should be a structure as returned by
%   loadBinaryPP, containing cell arrays `groups` and `thresholds`. If the
%   structure also contains a `subjects` field, the plot will display
%   different subjects in different colors. `pred` can be a vector or a
%   cell array of vectors giving the predicted thresholds with the same
%   ordering as that used in `meas`. When there are several prediction
%   vectors available, they will be plotted in different colors. The
%   predictions can also be absent (set `pred` to an empty matrix).
%
%   The function will ignore single-group entries.
%
%   Options:
%    'exclude'
%       Cell array of directions to exclude.
%    'predlabels'
%       Labels to show for the predictions.
%    'predellipses'
%       If true, draw ellipses for the predictions. Otherwise use a scatter
%       plot. If the predictions don't lie on an ellipse, a scatter plot is
%       used instead.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('exclude', {'A_1'}, @(c) iscell(c) && (isvector(c) || isempty(c)));
parser.addParameter('predlabels', {}, @(c) iscell(c) && (isvector(c) || isempty(c)));
parser.addParameter('predellipses', true, @(b) islogical(b) && isscalar(b));

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
edirs0 = {'A_1', 'AC_1_1', 'AB_1_1', 'AD_1_1', 'BC_1_1', ...
    'ABC_1_1_1', 'BCD_1_1_1', 'ABD_1_1_1', 'ACD_1_1_1', 'ABCD_1_1_1_1'};
keep_mask = ~ismember(edirs0, params.exclude);
edirs = edirs0(keep_mask);

% figure out which coordinate pairs we'll be looking at
edirs2 = cell(length(edirs), length(edirs));
for i = 1:length(edirs)
    for j = 1:i-1
        edirs2{i, j} = [edirs{i} ';' edirs{j}];
    end
end

% normalize group names to minimize the possibility for trivial mismatches
edirs2_norm = normalizegroup(edirs2);
meas_groups_norm = normalizegroup(meas.groups);

% read off sensitivities from predictions and measurements

% set up subjects
if isfield(meas, 'subjects')
    subjects = unique(meas.subjects);
    n_subjects = length(subjects);
else
    subjects = {'subj'};
    n_subjects = 1;
end

% figure out the colors
select_half = @(m) m(1:end/2, :);
% select_second_half = @(m) m(end/2+1:end, :);
% pred_colors = cool(size(length(pred), 1));
% pred_colors = select_second_half(cool(2*size(length(pred), 1)));
pred_colors = select_half(parula(2*size(length(pred), 1)));
if n_subjects > 1
    meas_colors = select_half(hot(2*n_subjects));
else
    meas_colors = [1 0 0];
end

% figure out all that we want to plot
to_plot = cell(size(edirs2));

for i = 1:size(edirs2, 1)
    for j = 1:i-1
        % find where to look in meas and pred
        crt_mask = strcmp(meas_groups_norm, edirs2_norm{i, j});
        
        % each element in to_plot is a structure
        crt_elem = repmat(struct, length(pred) + n_subjects, 1);
        
        % add the measurements
        for k = 1:n_subjects
            crt_subject = subjects{k};
            if isfield(meas, 'subjects')
                crt_submask = crt_mask & strcmp(meas.subjects, crt_subject);
            else
                crt_submask = crt_mask;
            end
            crt_dirs = meas.directions(crt_submask);
            crt_meas = meas.thresholds(crt_submask);
            crt_locs = binary2to1(binaryrec(crt_meas, crt_dirs));
            
            crt_elem(k).type = 'data';
            % j should map to x-position, i should map to y --> fliplr
            crt_elem(k).locs = fliplr(crt_locs);
            crt_elem(k).idx = k;
            crt_elem(k).color = meas_colors(k, :);
        end
        
        % add the predictions
        % use the last submask from the measurements -- this assumes that
        % all the predictions are the same for all subjects!
        for k = 1:length(pred)
            crt_pred = pred{k}(crt_submask);
            crt_locs = binary2to1(binaryrec(crt_pred, crt_dirs));
            
            crt_elem(n_subjects + k).type = 'pred';
            % j should map to x-position, i should map to y --> fliplr
            crt_elem(n_subjects + k).locs = fliplr(crt_locs);
            crt_elem(n_subjects + k).idx = k;
            crt_elem(n_subjects + k).color = pred_colors(k, :);
        end
        
        to_plot{i, j} = crt_elem;
    end
end

% fit ellipses, find major axes
to_plot_scales = zeros(size(to_plot));
for i = 1:size(edirs2, 1)
    for j = 1:i-1
        crt_elem = to_plot{i, j};
        for k = 1:length(crt_elem)
            crt_elem(k).ellipse = inv(fit_ellipse(crt_elem(k).locs));
            crt_elem(k).maj_axis = max(eig(crt_elem(k).ellipse));
        end
        
        to_plot_scales(i, j) = 0.25/max([crt_elem.maj_axis]);
        to_plot{i, j} = crt_elem;
    end
end

% finally, draw
was_hold = ishold;
hold on;

for i = 1:size(edirs2, 1)
    for j = 1:i-1
        crt_scaling = to_plot_scales(i, j);
        crt_elem = to_plot{i, j};
        for k = 1:length(crt_elem)
            crt_locs = crt_scaling*crt_elem(k).locs;
            crt_color = crt_elem(k).color;
            if strcmp(crt_elem(k).type, 'pred')
                done = false;
                if params.predellipses
                    done = ellipse(j, 10 - i, crt_scaling*crt_elem(k).ellipse, ...
                        'color', crt_color, 'linewidth', 1);
                end
                if ~done
                    scatter(j + crt_locs(:, 1), 10 - i + crt_locs(:, 2), [], crt_color, '.');
                end
            else
                ellipse(j, 10 - i, crt_scaling*crt_elem(k).ellipse, ...
                    'color', crt_color, 'linewidth', 0.5);
            end
        end
    end
end

% make the figure look good
axis equal;
xlim([0 1 + size(edirs2, 1)]);
ylim([0 1 + size(edirs2, 2)]);

% labels for the coordinates
c_labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% ticks
set(gca, 'xtick', (1:size(edirs2, 1)), 'xticklabel', c_labels(keep_mask));
set(gca, 'ytick', (1:size(edirs2, 2)), 'yticklabel', fliplr(c_labels(keep_mask)));

% title and legend
title('Precision matrix comparisons between natural images and psychophysics');

crt_elem = to_plot{2, 1};
legend_labels = cell(1, length(crt_elem));
for k = 1:length(crt_elem)
    if strcmp(crt_elem(k).type, 'pred')
        if isempty(params.predlabels)
            crt_label = ['pred ' int2str(crt_elem(k).idx)];
        else
            crt_label = params.predlabels(crt_elem(k).idx);
        end
    else
        if n_subjects > 1
            crt_label = ['data ' subjects{crt_elem(k).idx}];
        else
            crt_label = 'data';
        end
    end
    legend_labels{k} = crt_label;
end

legend(legend_labels);

% beautification
beautifygraph;

set(gca, 'xminortick', 'off', 'yminortick', 'off');

if ~was_hold
    hold off;
end

end