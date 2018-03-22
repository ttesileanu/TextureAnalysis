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

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('exclude', false(size(meas.thresholds)), @(b) islogical(b) && isvector(b) && length(b) == length(meas.thresholds));
parser.addParameter('errorbars', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('errcolor', [0.5 0.5 0.5], @(v) isvector(v) && isnumeric(v) && length(v) == 3);
parser.addParameter('ellipses', false, @(b) islogical(b) && isscalar(b));

if nargin == 1 && strcmp(pred, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% figure out whether we're doing errorbars
have_errors = params.errorbars & isfield(meas, 'threshold_intervals');

% figure out which groups we're showing
mask = ~params.exclude;
masked_pred = pred(mask);
masked_meas = meas.thresholds(mask);
if have_errors
    masked_meas_int = meas.threshold_intervals(mask, :);
end

masked_groups = meas.groups(mask);
masked_dirs = meas.directions(mask);

unique_groups = flipud(unique(masked_groups(:)));

% use the MatrixPlotter to make a figure containing a matrix of plots
plotter = MatrixPlotter(length(unique_groups));
while plotter.next
    i = plotter.index;

    hold on;
    
    % draw the background simplex triangle and unit and 1/2 radius circle,
    % for orientation
    drawTernaryTriangle;
    
    % locate measurements and predictions for the current group
    group_mask = strcmp(masked_groups, unique_groups{i});
    group_dirs = masked_dirs(group_mask);
    group_pred = masked_pred(group_mask);
    group_meas = masked_meas(group_mask);
    if have_errors
        group_meas_int = masked_meas_int(group_mask, :);
    end
    
    % convert predictions, thresholds, and threshold intervals to the 2d
    % coordinates necessary for plotting
    projected_pred = ternary3to2(ternaryrec(group_pred, group_dirs));
    projected_meas = ternary3to2(ternaryrec(group_meas, group_dirs));
    if have_errors
        projected_meas_lo = ternary3to2(ternaryrec(group_meas_int(:, 1), group_dirs));
        projected_meas_hi = ternary3to2(ternaryrec(group_meas_int(:, 2), group_dirs));
    end
    
    % draw error bars for measured data, if asked to
    if have_errors
        % don't plot errorbars that are not finite
        err_mask = all(isfinite(projected_meas_lo), 2) & all(isfinite(projected_meas_hi), 2);
        
        % each column is a plot
        plot([projected_meas_lo(err_mask, 1)' ; projected_meas_hi(err_mask, 1)'], ...
             [projected_meas_lo(err_mask, 2)' ; projected_meas_hi(err_mask, 2)'], ...
             'color', params.errcolor, 'linewidth', 0.5);
    end
    
    % draw measured data
    if have_errors
        % draw points without error bars as small open circles
        circle_mask = ~err_mask;
    else
        circle_mask = true(size(group_meas));
    end
    h_meas_cross = plot(projected_meas(~circle_mask, 1), projected_meas(~circle_mask, 2), ...
        'kx', 'linewidth', 1);
    h_meas_circ = plot(projected_meas(circle_mask, 1), projected_meas(circle_mask, 2), ...
        'ko', 'linewidth', 1, 'markersize', 3);
    
    % draw predictions
    h_pred = plot(projected_pred(:, 1), projected_pred(:, 2), ...
        'r.');
    
    if ~isempty(h_meas_cross)
        h_meas = h_meas_cross(1);
    else
        h_meas = h_meas_circ(1);
    end
    
    legend([h_meas h_pred(1)], 'measurements', 'predictions');
    
    % draw ellipses
    if params.ellipses
        % measurements
        crt_mask = all(isfinite(projected_meas), 2);
        if sum(crt_mask) > 2
            % remove non-finite entries
            meas_M = fit_ellipse(projected_meas(crt_mask, :));
            ellipse(0, 0, inv(meas_M), 'color', [0.5 0.5 0.5]);
        end
        
        % predictions
        crt_mask = all(isfinite(projected_pred), 2);
        if size(projected_pred, 1) > 2
            pred_M = fit_ellipse(projected_pred(crt_mask, :));
            ellipse(0, 0, inv(pred_M), 'color', [1 0.6 0.6]);
        end
    end
    
    % arrange plot sizes and labels
    axis equal;
    
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5]);
    
    xlim([-2 2]);
    ylim([-2 2]);
    
    title(unique_groups{i});
    
%    beautifygraph;
    set(gca, 'box', 'on', 'xminortick', 'on', 'yminortick', 'on', 'linewidth', 1);
end

preparegraph;

end