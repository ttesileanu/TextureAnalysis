function handles = makeMultiscalePCAPlot(pcs, pc_explained, N_values, varargin)
% makeMultiscalePCAPlot Plot principal components obtained from texture
% analysis at multiple scales.
%   makeMultiscalePCAPlot(pcs, pc_explained, N_values) plots the principal
%   components contained in the columns of the `pcs` matrix, splitting each
%   component into groups of 10 coefficients corresponding to the different
%   values of the block averaging factor `N` contained in `N_values`.
%   `N_values` and the vector of percent variance explained `pc_explained`
%   are just used for labeling the plots. The number of principal
%   components displayed is controlled by the 'maxPC' option (see below).
%
%   makeMultiscalePCAPlot({pc_scale_1, ..., pc_scale_K}, ...
%       {pc_explained_scale_1, ..., pc_explained_scale_K}, ...
%       N_values)
%   makes the plot when the principal component analysis was done
%   separately at each scale.
%
%   handles = makeMultiscalePCAPlot(...) returns a structure of handles ot
%   the various graphics elements drawn by the function.
%
%   makeMultiscalePCAPlot('defaults') display current default values for
%   the parameters.
%
%   Options:
%    'aspect'
%       Approximate aspect ratio for each bar plot.
%    'edges'
%       Pair [edge_left, edge_top] indicating size of edge left on the left
%       (for displaying values of `N`) and on the top (for displaying PC
%       index). Units are inches.
%    'maxPC'
%       Number of principal components to show.
%    'maxSize'
%       Pair [width, height] indicating maximum figure size (in inches).

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('aspect', 1.8, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x > 0);
parser.addParameter('edges', [0.5, 16/72], @(v) isnumeric(v) && isvector(v) && length(v) == 2 && ...
    isreal(v) && all(v > 0));
parser.addParameter('maxPC', 12, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x > 0);
parser.addParameter('maxSize', [19 11.5], @(v) isnumeric(v) && isvector(v) && length(v) == 2 && ...
    isreal(v) && all(v > 0));

% parse
parser.parse(varargin{:});
params = parser.Results;

% if asked, display defaults and exit
if ~iscell(pcs) && strcmp(pcs, 'defaults') && nargin == 1
    disp(params);
    return;
end

% handle situation in which PCA was done separately at each scale
if iscell(pcs)
    separate_pc = true;
    
    % just create a matrix `pcs` by concatenating individual PCs
    pcs = cell2mat(pcs(:));
else
    separate_pc = false;
end

% make sure maxPC isn't too large
if params.maxPC > size(pcs, 2)
    warning([mfilename ':badmaxPC'], ['Requested maxPC value (' int2str(params.maxPC) ...
        ') is larger than the number of available PCs (' int2str(size(pcs, 2)) ...
        '). Setting maxPC to ' int2str(size(pcs, 2)) '.']);
    params.maxPC = size(pcs, 2);
end

% figure out figure size
handles.fig = figure;
handles.fig.Units = 'inches';

proposed_y_each = params.maxSize(2)/length(N_values);
proposed_x_from_y = proposed_y_each*params.aspect*params.maxPC;
if proposed_x_from_y > params.maxSize(1)
    fig_x = params.maxSize(1);
    fig_y = params.maxSize(1)*length(N_values)/params.maxPC/params.aspect;
else
    fig_y = params.maxSize(2);
    fig_x = params.maxSize(2)*params.maxPC*params.aspect/length(N_values);
end
handles.fig.Position = [0.5 0.5 fig_x fig_y];

% texture statistics labels
labels = {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'};

% figure out figure edges in normalized coordinates
edge_left_normalized = params.edges(1) / handles.fig.Position(3);
edge_top_normalized = params.edges(2) / handles.fig.Position(4);

% draw
for j = 1:length(N_values)
    % the indices within the PC corresponding to the current scalr
    idxs = (1 + (j-1)*10):j*10;
    for i = 1:params.maxPC
        % figure out axis placement on the figure
        ax = axes;
        handles.barplots{i, j} = ax;
        
        crt_width = (1 - edge_left_normalized)/params.maxPC;
        crt_left_edge = edge_left_normalized + crt_width*(i-1);
        crt_height = (1 - edge_top_normalized)/length(N_values);
        crt_bottom_edge = crt_height*(length(N_values) - j);
        set(ax, 'units', 'normalized', 'outerposition', ...
            [crt_left_edge, crt_bottom_edge, ...
             crt_width, crt_height]);
        
        % normalize principal component so that largest (in absolute value)
        % entry has magnitude 1
        crt_pc = pcs(idxs, i);
        crt_pc_norm = crt_pc / max(abs(crt_pc) + eps);
        
        % draw and set limits in a consistent way
        bar(crt_pc_norm, 'edgecolor', 'none');
        xlim([0.5 10.5]);
        ylim([-1 1]);
        
        % use the fancy labels for the tick marks
        set(ax, 'xtick', 1:10, 'xticklabel', labels);
        
        beautifygraph;
        
        if separate_pc && iscell(pc_explained)
            % if the principal components are given separately, write
            % percent explained within each plot
            xl = xlim;
            yl = ylim;
            text(xl(2), yl(2), [int2str(round(pc_explained{j}(i))) '%'], ...
                'horizontalalignment', 'right', ...
                'verticalalignment', 'cap');
        end

        % for the first row of plots, print out a label identifying the
        % rank of the pricipal component, and the amount of variance
        % explained (unless we have variance separately for each PC)
        if j == 1
            ax_top = axes;
            handles.top_labels{i} = ax_top;
            
            set(ax_top, 'units', 'normalized', 'position', ...
                [crt_left_edge, 1 - edge_top_normalized, ...
                crt_width, edge_top_normalized]);
            axis off;
            if ~separate_pc || ~iscell(pc_explained)
                text(0.5, 0.5, ['PC' int2str(i) '  ' int2str(round(pc_explained(i))) '%']);
            else
                text(0.5, 0.5, ['PC' int2str(i)]);
            end
        end
    end
    
    % next to the first column of plots, display the corresponding block
    % averaging factor (N)
    ax_left = axes;
    handles.left_labels{j} = ax_left;
    
    set(ax_left, 'units', 'normalized', 'position', ...
        [0, crt_bottom_edge, edge_left_normalized, crt_height]);
    axis off;
    text(0.5, 0.6, ['N=' int2str(N_values(j))], 'horizontalalignment', 'center');
end

% set figure position and bounding box consistently so that saving to PDF
% works well
preparegraph;

end