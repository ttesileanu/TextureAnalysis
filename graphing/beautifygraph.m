function beautifygraph(varargin)
% BEAUTIFYGRAPH Improve the quality of a graph.
%   BEAUTIFYGRAPH improves the quality of the current active axes by
%   changing font styles and sizes and tick directions.
%
%   BEAUTIFYGRAPH(haxes) is applied to the axes identified by the given
%   handle.
%
%   Options:
%    'changefonts' <b>
%       True to change fonts.
%    'fontscale' <x>
%       A scale factor for all the font sizes.
%    'tickdir'
%    'box'
%    'linewidth'
%       Options to set for the axes.
%    'ticks'
%    'minorticks'
%       Set to 'on' or 'off' to either show or not show the (minor) ticks.
%    'fontnames'
%       Set fonts to use. This should be a pair, {axes font, title/labels
%       font}.
%    'titleweight'
%       Set to 'bold' to display the title in bold letters, set to 'normal'
%       otherwise.
%    'titlesize'
%       Font size for title. This is further scaled by 'fontscale'.
%    'labelsize'
%       Font size for labels. This is further scaled by 'fontscale'.
%    'ticksize'
%       Font size for tick labels. This is further scaled by 'fontscale'.
%    'ticklabels'
%       Set to false to remove all tick labels.
%    'noaxes'
%       Set to true to not display the axes. This includes the background.

% This function is based on the suggestions given at
% blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs

% parse the optional arguments that scattefit cares about
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('axes', []);
parser.addParameter('fontscale', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('changefonts', false, @(b) isscalar(b) && islogical(b));
parser.addParameter('tickdir', 'out', @(s) ismember(s, {'in', 'out'}));
parser.addParameter('box', 'off', @(s) ismember(s, {'on', 'off'}));
parser.addParameter('linewidth', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('minorticks', 'on', @(s) ismember(s, {'on', 'off'}));
parser.addParameter('ticks', 'on', @(s) ismember(s, {'on', 'off'}));
parser.addParameter('fontnames', {'Helvetica', 'AvantGarde'}, ...
    @(c) iscell(c) && all(cellfun(@(s) ischar(s) && isvector(s), c)));
parser.addParameter('titleweight', 'bold', @(s) ismember(s, {'bold', 'normal'}));
parser.addParameter('titlesize', 14, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('labelsize', 12, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('ticksize', [], @(x) isscalar(x) && isnumeric(x));
parser.addParameter('ticklabels', true, @(b) isscalar(b) && islogical(b));
parser.addParameter('noaxes', false, @(b) isscalar(b) && islogical(b));

% show defaults
if length(varargin) == 1 && strcmp(varargin{1}, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% use current axes if none provided
if isempty(params.axes)
    params.axes = gca;
end

% find handles for relevant elements
htitle = get(params.axes, 'title');
hxlabel = get(params.axes, 'xlabel');
hylabel = get(params.axes, 'ylabel');

% change fonts if necessary
if params.changefonts
    set(params.axes, 'fontname', params.fontname{1});
    set([htitle, hxlabel, hylabel], 'fontname', params.fontname{2});
end

% set size of tick labels
if ~isempty(params.ticksize)
    xax = get(params.axes, 'xaxis');
    set(xax, 'fontsize', params.ticksize*params.fontscale);
    
    yax = get(params.axes, 'yaxis');
    set(yax, 'fontsize', params.ticksize*params.fontscale);
end

% set sizes of axis labels and title
set([hxlabel, hylabel], 'fontsize', params.labelsize*params.fontscale);
set(htitle, 'fontsize', params.titlesize*params.fontscale, 'fontweight', params.titleweight);

% update some other options
set(params.axes, ...
    'box',          params.box, ...
    'tickdir',      params.tickdir, ...
    'xminortick',   params.minorticks, ...
    'yminortick',   params.minorticks, ...
    'linewidth',    params.linewidth);
% disable ticks if asked
if strcmp(params.ticks, 'off')
    set(params.axes, 'ticklength', [0 0]);
end
% disable tick labels if asked
if ~params.ticklabels
    set(params.axes, 'xtick', [], 'ytick', []);
end
% remove the axes if asked
if params.noaxes
    axis off;
end

end
