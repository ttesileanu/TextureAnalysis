function ternaryPredictionMatchPerGroup(pred, meas, varargin)
% ternaryPredictionMatchPerGroup Make plots for each texture group
% comparing predicted to measured thresholds.
%   ternaryPredictionMatchPerGroup(pred, meas) makes a plot for each
%   texture group present in the measurement structure `meas`. This
%   structure should contain `groups`, `directions`, and `thresholds`. If
%   `threshold_intervals` are present, they are also used (unless error
%   bars are disabled; see options below). The predicted thresholds `pred`
%   should be given as a vector with the same ordering as that in the
%   measurements structure.
%
%   Options:
%    'exclude'
%       Boolean mask showing which measurements to exclude from the plots.
%    'errorbars'
%       Set to true to draw error bars, false otherwise.
%    'errcolor'
%       3-component RGB vector giving the color of the error bars.
%    'ellipses'
%       Set to true to show best-fit ellipses, false otherwise.
%    'colorfct'
%       Colormap function giving the colors to use when there are several
%       subjects. colorfct(n) should return an n x 3 matrix of RGB colors.
%       Set to empty to disable per-subject coloring.
%    'multi'
%       How to deal with mixed groups. If set to `false`, mixed groups are
%       ignored. If set to `true`, all groups are drawn. It can also be a
%       vector giving the multiplicities to draw, for instance [1] shows
%       only single groups, [2] shows pairs, and [1, 2] shows both. Note
%       that for now only single groups and pairs are supported.
%    'beautifyopts'
%       Options to pass to beautifygraph.
%    'triangleopts'
%       Options to pass to drawTernaryTriangle.
%    'limits'
%       The largest coordinate value to show in absolute value and in
%       either direction. This can be a pair of numbers, in which case the
%       first applies to single groups, and the second to mixed groups.
%    'plotter_opts'
%       Options to pass to MatrixPlotter.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

if nargin == 1 && strcmp(pred, 'defaults')
    show_defaults = true;
    meas.thresholds = [];
else
    show_defaults = false;
end

parser.addParameter('exclude', false(size(meas.thresholds)), @(b) islogical(b) && isvector(b) && length(b) == length(meas.thresholds));
parser.addParameter('errorbars', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('errcolor', [0.5 0.5 0.5], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('ellipses', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('colorfct', @parula, @(f) isempty(f) || isa(f, 'function_handle'));
parser.addParameter('multi', false, @(b) (islogical(b) && isscalar(b)) | ...
    (isnumeric(b) && isvector(b)));
parser.addParameter('beautifyopts', {'box', 'on', 'tickdir', 'in', 'titlesize', 12}, ...
    @(c) iscell(c) && (isempty(c) || isvector(c)));
parser.addParameter('limits', [2 1], @(v) isnumeric(v) && isvector(v) && ismember(length(v), [1 2]));
parser.addParameter('triangleopts', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('plotter_opts', {}, @(c) iscell(c) && isvector(c));

if show_defaults
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% figure out what to do with multiple groups
if islogical(params.multi)
    if params.multi
        params.multi = [1, 2];
    else
        params.multi = 1;
    end
end

% handle 1-element limits
if length(params.limits) == 1
    params.limits = repmat(params.limits, 1, 2);
end

% figure out whether we're doing errorbars
have_errors = params.errorbars & isfield(meas, 'threshold_intervals');

% very basic input cleaning: remove whitespace from group names
meas.groups = cellfun(@(s) s(~isspace(s)), meas.groups, 'uniform', false);

% figure out which entries are tuples, and of how many elements
multiplicities = 1 + cellfun(@(s) sum(s == ';'), meas.groups);
% figure out what to keep
multi_mask = ismember(multiplicities, params.multi);

% figure out which groups we're showing
mask = ~params.exclude & multi_mask;
masked_pred = pred(mask);
masked_meas = meas.thresholds(mask);
if have_errors
    masked_meas_int = meas.threshold_intervals(mask, :);
end

masked_groups = meas.groups(mask);
masked_dirs = meas.directions(mask);

if ~isempty(params.colorfct) && isfield(meas, 'subjects')
    masked_subj = meas.subjects(mask);
    unique_subj = unique(masked_subj);
    n_subjects = length(unique_subj);
    have_colors = n_subjects > 1;
    subj_colors = params.colorfct(n_subjects);
else
    have_colors = false;
    n_subjects = 1;
end

unique_groups = sortgroups(unique(masked_groups(:)));

% use the MatrixPlotter to make a figure containing a matrix of plots
plotter = MatrixPlotter(length(unique_groups), params.plotter_opts{:});
plotter.ax_aspect = 1;
while plotter.next
    i = plotter.index;

    hold on;

    % figure out the group(s) we're working with
    crt_tuple = unique_groups{i};
    crt_groups = strsplit(crt_tuple, ';');
    crt_multi = length(crt_groups);
    
    if crt_multi == 1
        % draw the background simplex triangle and unit and 1/2 radius circle,
        % for orientation
        drawTernaryTriangle(params.triangleopts{:});
    elseif crt_multi == 2
        % draw axes through origin and guiding circles
        drawTernaryMixedBackground(crt_groups{:});
    else
        error([mfilename ':badmulti'], 'Only 1- or 2-tuples of groups are supported.');
    end
    
    % locate measurements and predictions for the current group
    group_mask = strcmp(masked_groups, crt_tuple);
    group_dirs = masked_dirs(group_mask);
    group_pred = masked_pred(group_mask);
    group_meas = masked_meas(group_mask);
    if have_errors
        group_meas_int = masked_meas_int(group_mask, :);
    end
    if have_colors
        group_subj = masked_subj(group_mask);
    end

    % convert predictions, thresholds, and threshold intervals to the 2d
    % coordinates necessary for plotting
    if crt_multi == 1
        projected_pred = ternary3to2(ternaryrec(group_pred, group_dirs));
        projected_meas = ternary3to2(ternaryrec(group_meas, group_dirs));
        if have_errors
            projected_meas_lo = ternary3to2(ternaryrec(group_meas_int(:, 1), group_dirs));
            projected_meas_hi = ternary3to2(ternaryrec(group_meas_int(:, 2), group_dirs));
        end
    elseif crt_multi == 2
        projected_pred = ternary6tomix2(ternaryrec(group_pred, group_dirs));
        projected_meas = ternary6tomix2(ternaryrec(group_meas, group_dirs));
        if have_errors
            projected_meas_lo = ternary6tomix2(ternaryrec(group_meas_int(:, 1), group_dirs));
            projected_meas_hi = ternary6tomix2(ternaryrec(group_meas_int(:, 2), group_dirs));
        end
    end
    
    % draw error bars for measured data, if asked to
    if have_errors
        % don't plot errorbars that are not finite
        err_mask = all(isfinite(projected_meas_lo), 2) & all(isfinite(projected_meas_hi), 2);

        if ~have_colors
            % each column is a plot
            plot([projected_meas_lo(err_mask, 1)' ; projected_meas_hi(err_mask, 1)'], ...
                 [projected_meas_lo(err_mask, 2)' ; projected_meas_hi(err_mask, 2)'], ...
                 'color', params.errcolor, 'linewidth', 0.5);
        else
             % figure out per-subject coloring
            err_meas_lo = projected_meas_lo(err_mask, :);
            err_meas_hi = projected_meas_hi(err_mask, :);
            err_subj = group_subj(err_mask);
            
            for k = 1:size(err_meas_lo, 1)
                subj_id = find(strcmp(unique_subj, err_subj{k}), 1);
                subj_color = subj_colors(subj_id, :);
                plot([err_meas_lo(k, 1) err_meas_hi(k, 1)], ...
                     [err_meas_lo(k, 2) err_meas_hi(k, 2)], ...
                     'color', mixcolor(subj_color, params.errcolor), ...
                     'linewidth', 0.5);
            end
         end
    end
    
    % draw measured data
    if have_errors
        % draw points without error bars as small open circles
        circle_mask = ~err_mask;
    else
        circle_mask = true(size(group_meas));
    end
    
    h_meas_cross = [];
    h_meas_circ = [];
    for k = 1:n_subjects
        if have_colors
            crt_subj = unique_subj(k);
            subj_mask = strcmp(group_subj, crt_subj);
            if sum(subj_mask) == 0
                continue;
            end
        
            crt_meas_color = subj_colors(k, :);
        else
            subj_mask = true(size(circle_mask));
            crt_meas_color = [0 0 0];
        end
        
        h_meas_cross0 = plot(projected_meas(~circle_mask & subj_mask, 1), ...
            projected_meas(~circle_mask & subj_mask, 2), ...
            'x', 'linewidth', 1, 'color', crt_meas_color);
        h_meas_circ0 = plot(projected_meas(circle_mask & subj_mask, 1), ...
            projected_meas(circle_mask & subj_mask, 2), ...
            'o', 'linewidth', 1, 'markersize', 3, 'color', crt_meas_color);
        if isempty(h_meas_cross)
            h_meas_cross = h_meas_cross0;
        end
        if isempty(h_meas_circ)
            h_meas_circ = h_meas_circ0;
        end
                    
        % draw ellipses
        if params.ellipses
            % measurements
            crt_mask = all(isfinite(projected_meas), 2);
            if sum(crt_mask) > 2
                % remove non-finite entries
                meas_M = fit_ellipse(projected_meas(crt_mask & subj_mask, :));
                ellipse(0, 0, inv(meas_M), 'color', mixcolor(crt_meas_color, [0.5 0.5 0.5]));
            end            
        end
    end
    
    if ~isempty(h_meas_cross)
        h_meas = h_meas_cross(1); %#ok<NASGU>
    else
        h_meas = h_meas_circ(1); %#ok<NASGU>
    end

    % draw predictions
    h_pred = plot(projected_pred(:, 1), projected_pred(:, 2), ...
        'r.'); %#ok<NASGU>
    
    % draw ellipses
    if params.ellipses
        % predictions
        crt_mask = all(isfinite(projected_pred), 2);
        if size(projected_pred, 1) > 2
            pred_M = fit_ellipse(projected_pred(crt_mask, :));
            ellipse(0, 0, inv(pred_M), 'color', [1 0.6 0.6]);
        end
    end
    
%     legend([h_meas h_pred(1)], 'measurements', 'predictions');
    
    % arrange plot sizes and labels
    axis equal;
    
    if crt_multi == 1
        xlim([-params.limits(1) params.limits(1)]);
        ylim([-params.limits(1) params.limits(1)]);
    else
        xlim([-params.limits(2) params.limits(2)]);
        ylim([-params.limits(2) params.limits(2)]);
    end
    
    title(crt_tuple);
    
    beautifygraph(params.beautifyopts{:});
end

fig = plotter.get_figure;
fig.Name = [fig.Name ' (crosses, open circles = PP, dots = NI)'];

preparegraph;

end

function color = mixcolor(color1, color2)
% Mix the two colors by adding them.

color = max(min(color1 + color2, 1), 0);

end