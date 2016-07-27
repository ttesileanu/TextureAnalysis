function [evProj, evFProj, evBProj, pcaData] = pcaForeVsBack(data, varargin)
% pcaForeVsBack Project foreground and background statistics onto principal
% components calculated from all patches.
%   [evProj, evFProj, evBProj] = pcaForeVsBack(data) finds the top 2
%   principal components of the texture statistics from data.ev, and
%   projects data.evF and data.evB on these components. The function then
%   shows the patches on a plot (unless 'plot' is false; see options below).
%
%   [evProj, evFProj, evBProj, pcaData] = pcaForeVsBack(data) also returns
%   a structure containing details about the PCA. pcaData.U, pcaData.S,
%   pcaData.V are the U, S, V matrices returned from the svd command
%   ('economy' size); while pcaData.mu is the mean of the patch statistics
%   from data.ev.
%
%   Options:
%    'alpha': float
%       Amount of transparency to use for the scatter plots.
%    'color': {colFore, colBack, colNeither}
%       Sets the colors to use for foreground, background, and patches that
%       are neither. Set any one of these to an empty string to supress
%       drawing that component.
%    'plot': bool
%       Set to false to suppress the plot.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('alpha', 0.5, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('colors', {'b', 'r', [0.8, 0.8, 0.8]}, @(c) iscell(c) && isvector(c));
parser.addParameter('plot', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% mean-center full patch results, find SVD
mu = mean(data.ev, 1);
evCentered = bsxfun(@minus, data.ev, mu);
[svdU, svdS, svdV] = svd(evCentered, 0);

evProj = evCentered*svdV;
evFProj = bsxfun(@minus, data.evF, mu)*svdV;
evBProj = bsxfun(@minus, data.evB, mu)*svdV;

pcaData.mu = mu;
pcaData.U = svdU;
pcaData.S = svdS;
pcaData.V = svdV;

if params.plot
    oldHold = ishold;
    hold on;
    
    plotOpts = {'filled'};
    
    if ~isempty(params.colors{3})
        scatter(evProj(:, 1), evProj(:, 2), [], params.colors{3}, plotOpts{:});
        h3 = plot(evProj(1, 1), evProj(1, 2), '.', 'color', params.colors{3}, 'markersize', 20);
        set(h3, 'visible', 'off');
    else
        h3 = [];
    end
    if ~isempty(params.colors{1})
        scatter(evFProj(:, 1), evFProj(:, 2), [], params.colors{1}, plotOpts{:});
        h1 = plot(evFProj(1, 1), evFProj(1, 2), '.', 'color', params.colors{1}, 'markersize', 20);
        set(h1, 'visible', 'off');
    else
        h1 = [];
    end
    if ~isempty(params.colors{2})
        scatter(evBProj(:, 1), evBProj(:, 2), [], params.colors{2}, plotOpts{:});
        h2 = plot(evBProj(1, 1), evBProj(1, 2), '.', 'color', params.colors{2}, 'markersize', 20);
        set(h2, 'visible', 'off');
    else
        h2 = [];
    end
    
    handles = [];
    if ~isempty(h1)
        handles = [handles h1];
    end
    if ~isempty(h2)
        handles = [handles h2];
    end
    if ~isempty(h3)
        handles = [handles h3];
    end
    
    alpha(params.alpha);
    xlabel('PC #1');
    ylabel('PC #2');
    
    legend(handles, {'foreground', 'background', 'unassigned'});
    
    beautifygraph;
    preparegraph;
    
    if ~oldHold
        hold off;
    end
end

end