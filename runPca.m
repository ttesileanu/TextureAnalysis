function [evProj, evObjProj, pcaData] = runPca(res, varargin)
% runPca Project statistics onto principal components.
%   [evProj, evObjProj] = pcaForeVsBack(res) finds the top 2 principal
%   components of the texture statistics from `res.ev`, and projects the
%   data onto these components. `evObjProj` is a `Map` object giving the
%   projections for patches restricted to each object id. The function
%   shows the patches on a plot (unless 'plot' is false; see options below).
%
%   [..., pcaData] = pcaForeVsBack(res) also returns a structure containing
%   details about the PCA. `pcaData.U`, `pcaData.S`, `pcaData.V` are the U,
%   S, V matrices returned from the `svd` command ('economy' size); while
%   `pcaData.mu` is the mean of the patch statistics from `res.ev`.
%
%   Options:
%    'alpha': float
%       Amount of transparency to use for the scatter plots.
%    'plot': bool
%       Set to false to suppress the plot.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('alpha', [], @(x) isnumeric(x) && isscalar(x));
parser.addParameter('plot', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% mean-center full patch results, find SVD
mu = mean(res.ev, 1);
evCentered = bsxfun(@minus, res.ev, mu);
[svdU, svdS, svdV] = svd(evCentered, 0);

evProj = evCentered*svdV;

allObjIds = unique(res.objIds);
evObjProj = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:length(allObjIds)
    objId = allObjIds(i);
    evObjProj(i) = evProj(res.objIds == objId, :);
end

pcaData.mu = mu;
pcaData.U = svdU;
pcaData.S = svdS;
pcaData.V = svdV;

if params.plot
    oldHold = ishold;
    hold on;
    
    plotOpts = {'filled'};
    if ~isempty(params.alpha)
        plotOpts = [plotOpts {'alpha', params.alpha}];
    end
    
    smartscatter(evProj(:, 1), evProj(:, 2), [], res.objIds, plotOpts{:});
    colormap('parula');
    xlabel('PC #1');
    ylabel('PC #2');
    
    handles = [];
    colors = parula(length(allObjIds));
    for i = 1:length(allObjIds)
        handles(i) = plot(evProj(1, 1), evProj(1, 2), '.', 'color', colors(i, :), 'markersize', 20); %#ok<AGROW>
        set(handles(i), 'visible', 'off');
    end
    
    legend(handles, arrayfun(@(i) ['Object ' int2str(i)], allObjIds, 'uniform', false));
    
    beautifygraph;
    preparegraph;
    
    if ~oldHold
        hold off;
    end
end

end