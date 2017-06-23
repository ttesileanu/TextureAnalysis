function plotStatisticsBars(evData, varargin)
% plotStatisticsBars Make bar plot showing the texture statistics.
%   plotStatisticsBars(evData) makes a bar plot showing the means and
%   standard errors for the given texture statistics.
%
%   plotStatisticsBars({evData1, evData2, ...}) plots several different
%   sets of statistics on the same axes.
%
%   plotStatisticsBars(evStruct), when evStruct is a structure that has
%   fields 'evF' and 'evB', is equivalent to
%   plotStatisticsBars({evStruct.evF, evStruct.evB}).
%
%   Options:
%    'aspect': float
%       Aspect ratio (ratio between horizontal and vertical sizes of the
%       plot).
%    'colors': cell array
%       The function cycles through these colors to draw the bars. The
%       colors can be given as Matlab strings, or RGB triplets.
%    'spacing': float
%       This sets the amount of space between the columns corresponding to
%       different statistics (as a fraction of the total column width).
%    'datasetSpacing': float
%       This sets the amount of space between the bars corresponding to
%       different data sets (as a fraction of the total column width).
%    'legend': cell array of string
%       Gives the names to use for the legend. Set to empty to have no
%       legend.
%    'type': str
%       This can be 'bar' or 'whisker'.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('aspect', 7/3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('colors', [], @(c) isempty(c) || (iscell(c) && isvector(c)));
parser.addParameter('spacing', 0.15, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('datasetSpacing', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('type', 'whisker', @(s) ischar(s) && isvector(s));
parser.addParameter('legend', {'foreground', 'background'}, @(c) isempty(c) || ...
    (iscell(c) && isvector(c)));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~iscell(evData)
    if isstruct(evData)
        evData = {evData.evF, evData.evB};
    else
        evData = {evData};
    end
end

if isempty(params.colors)
    if strcmp(params.type, 'whisker')
        params.colors = {'k', 'r', 'b'};
    else
        params.colors = {'k', 'w', 'r', 'b'};
    end
end

nEv = numel(evData);
if strcmp(params.type, 'whisker')
    barWidth = 0.15;
else
    barWidth = (1 - params.spacing - (nEv - 1)*params.datasetSpacing)/nEv;
end
datasetSkip = barWidth + params.datasetSpacing;

figure;
hold on;
handles = [];
for i = 1:nEv
    skip = (i - 1)*datasetSkip;
    iCycle = 1 + mod(i-1, length(params.colors));
    
    if strcmp(params.type, 'whisker')
        hAll = boxplot(evData{i}, 'positions', 1+skip:10+skip, ...
            'colors', params.colors{iCycle}, 'symbol', '', 'boxstyle', 'filled', ...
            'medianstyle', 'target');
        h = findall(hAll, 'tag', 'Box');
        handles = [handles; h(1)]; %#ok<AGROW>
    elseif strcmp(params.type, 'bar')
        data = mean(evData{i}, 1);
        bar(1+skip:10+skip, data, barWidth, 'facecolor', params.colors{iCycle});
        if size(evData{i}, 1) > 1
            err = std(evData{i}, 1);
            
            errorbar(1+skip:10+skip, data, err, 'color', [0.5, 0.5, 0.5], 'linestyle', 'none');
        end
    else
        error([mfilename ':badtype'], 'Unrecognized plot type.');
    end
end

if strcmp(params.type, 'whisker')
    shift = 0.1;
else
    shift = 0.5*(1 - barWidth);
end
set(gca, 'xtick', 1+shift:10+shift, 'xticklabel', ...
    {'\gamma', '\beta_{|}', '\beta_{--}', '\beta_{\\}', '\beta_{/}', ...
    '\theta_{\lceil}', '\theta_{\rfloor}', '\theta_{\rceil}', '\theta_{\lfloor}', '\alpha'});

% set aspect ratio
oldPos = get(gcf, 'position');
set(gcf, 'position', [oldPos(1:2) oldPos(3) oldPos(3)/params.aspect]);

if strcmp(params.type, 'whisker')
    xlim([0.5, 11]);
    set(gca, 'xaxislocation', 'origin', 'box', 'off');
    
    % a hack to take the x-axis labels off the axis
    delta = 0.01;
    Xt = get(gca, 'xtick');
    Yl = ylim;
    t = text(Xt, Yl(1)*ones(1, length(Xt))-delta, get(gca, 'xticklabel'));
    set(t, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    set(gca,'XTickLabel','');
    
    for k = 1:length(t)
        set(t(k), 'units', 'normalized');
    end
    
    if ~isempty(params.legend)
        legend(handles, params.legend);
    end
else
    if ~isempty(params.legend)
        legend(params.legend);
    end
end

beautifygraph;
preparegraph;

set(gca, 'xminortick', 'off');

end