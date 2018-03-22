function showPredictionMatch(pred, meas, groups, varargin)
% showPredictionMatch Make a scatter plot comparing predicted and observed
% values, grouped by color.
%   showPredictionMatch(pred, meas, groups) makes a scatter plot showing
%   how the measured values `meas` compare to predicted values `pred`,
%   grouping them by color based on the string values in the cell array
%   `groups`. Entries that are either infinite or NaN in either `pred` or
%   `meas` are ignored.
%
%   Options:
%    'exclude'
%       Boolean mask of entries to exclude.
%    'measintervals'
%       Error intervals for the measurements, as an n x 2 matrix, where n
%       is the number of measurements.
%    'subjects'
%       Cell array identifying the subject for which each measurement was
%       obtained. Used when 'subjavg' is true.
%    'subjavg'
%       Show only per-subject average for each group.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('exclude', false(size(meas)), @(b) islogical(b) && isvector(b) && length(b) == length(meas));
parser.addParameter('measintervals', [], @(m) isnumeric(m) && size(m, 1) == length(meas) && size(m, 2) == 2);
parser.addParameter('subjavg', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('subjects', {}, @(c) iscell(c) && length(c) == length(meas));

if nargin == 1 && strcmp(pred, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% exclude non-finite elements, and those explicitly excluded by the user
mask = isfinite(pred) & isfinite(meas) & ~params.exclude;

masked_groups = groups(mask);
masked_pred = pred(mask);
masked_meas = meas(mask);

% we will need hold on, but make sure to return it to whatever it was
washold = ishold;
hold on;

% figure out what colors to use -- need to know how many groups we have
unique_groups = fliplr(unique(masked_groups));
n_unique_groups = length(unique_groups);
group_colors = parula(n_unique_groups);

axis equal;

if ~params.subjavg
    colors = zeros(sum(mask), 3);
    for i = 1:n_unique_groups
        for k = 1:3
            colors(strcmp(masked_groups, unique_groups{i}), k) = group_colors(i, k);
        end
    end
    
    % draw the error bars first
    err_pos = params.measintervals(:, 2) - meas;
    err_neg = meas - params.measintervals(:, 1);
    
    h = errorbar(masked_pred, masked_meas, err_neg(mask), err_pos(mask), ...
        'marker', 'none', 'color', [0.5 0.5 0.5], 'linestyle', 'none');
    h.CapSize = 0;
    
    % draw colored lines connecting al the dots within each group
    for i = 1:n_unique_groups
        sub_mask_group = strcmp(masked_groups, unique_groups{i});
        sub_pred = masked_pred(sub_mask_group);
        sub_meas = masked_meas(sub_mask_group);
        
        crt_color = 0.5 + 0.5*group_colors(i, :);
        
        for j = 1:length(sub_pred)
            for k = j+1:length(sub_pred)
                line(sub_pred([j k]), sub_meas([j k]), 'linewidth', 0.5, ...
                    'color', crt_color);
            end
        end
    end
    
    % make the scatter plot of measurements vs. predictions
    smartscatter(masked_pred, masked_meas, 'color', colors, 'density', false);
    
    all_data = [masked_meas(:) ; masked_pred(:)];    
else
    % check that we have subjects names
    if isempty(params.subjects)
        params.subjects = repmat({'default'}, length(meas), 1);
    end
    
    masked_subjects = params.subjects(mask);
    
    % take averages by group, by subject
    groupSubjPairs = {};
    avg_pred = [];
    avg_meas = [];
    colors = [];
    for i = 1:n_unique_groups
        crt_group = unique_groups{i};
        group_mask = strcmp(masked_groups, crt_group);
        crt_subjects = masked_subjects(group_mask);
        for j = 1:length(crt_subjects)
            sub_mask = group_mask & strcmp(masked_subjects, crt_subjects{j});
            crt_avg_pred = mean(masked_pred(sub_mask));
            crt_avg_meas = mean(masked_meas(sub_mask));
            avg_pred = [avg_pred crt_avg_pred]; %#ok<AGROW>
            avg_meas = [avg_meas crt_avg_meas]; %#ok<AGROW>
            groupSubjPairs = [groupSubjPairs {{crt_group crt_subjects{j}}}]; %#ok<AGROW>
            colors = [colors ; group_colors(i, :)]; %#ok<AGROW>
        end
    end
    
    % now draw
    smartscatter(avg_pred, avg_meas, 'density', false, 'color', colors);
    
    all_data = [avg_meas(:) ; avg_pred(:)];    
end

% use equal scaling on the x- and y-axes, and figure out an appropriate range
t_min = min(all_data);
t_max = max(all_data);

ylim([0.9*t_min 1.1*t_max]);

% draw the identity line, together with the correlation coefficient between
% (average) measurements and (average) predictions
if ~params.subjavg
    drawfitline(masked_pred, masked_meas, 'line', [1 0], ...
        'style', {'--k'}, 'legendloc', 'northwest');
else
    drawfitline(avg_pred, avg_meas, 'line', [1 0], ...
        'style', {'--k'}, 'legendloc', 'northwest');
end

% set labels
xlabel('Predicted thresholds');
ylabel('Actual thresholds');

if ~params.subjavg
    % label the dots using their group names
    for i = 1:length(meas)
        if ~mask(i)
            continue;
        end
        
        text(pred(i)+0.01, meas(i), groups{i}, 'fontsize', 6);
    end
end

beautifygraph;

if ~washold
    hold off;
end

end